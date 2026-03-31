# =============================================================================
# Hawkes M1 - 0 Change Points
# Baseline: Weibull
# Kernel:   exponential, eta and delta global
# Model:    lambda(t) = (alpha/beta)*(t/beta)^(alpha-1) + eta * A(t)
#           eta = n * delta,  n ~ Unif(0,1),  delta ~ Unif(0,3)
# =============================================================================

library(xtable)
library(rjags)
library(lattice)
library(stats)
library(bayesplot)
library(loo)
library(R2jags)

options(mc.cores = 3)
set.seed(42)

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
t <- 1096  # Total days

di <- list(excePM25col_Bog, excePM25col_Med,
           excePM25who_Bog, excePM25who_Med)

d_names <- c("PM25col_Bog", "PM25col_Med",
             "PM25who_Bog", "PM25who_Med")

# -----------------------------------------------------------------------------
# MCMC settings
# Slightly more iterations than NHPP-0CP due to additional parameters n, delta
# -----------------------------------------------------------------------------
niter   <- 100000
nburnin <- 20000
nthin   <- 10

# Constant C for zeros trick: must satisfy phi[k] = -log(lambda_k) + l2 + l3 + C >= 0
# C = 10 is conservative for typical exceedance intensities
C <- 10

# -----------------------------------------------------------------------------
# Estimation loop
# Uncomment the line below before first run:
samples_hawkes0 <- list()
# -----------------------------------------------------------------------------

t1 <- Sys.time()

for (i in 1:4) {

  d <- di[[i]]
  K <- length(d)

  line_data <- list(
    "K" = K,
    "T" = t,
    "d" = d,
    "C" = C
  )

  samples_hawkes0[[i]] <- jags(
    data               = line_data,
    parameters.to.save = c("alpha", "beta", "n", "delta", "eta", "A", "phi2"),
    n.chains           = 2,
    n.iter             = niter,
    n.burnin           = nburnin,
    n.thin             = nthin,
    model.file         = "02_modelos_JAGS/bug_files/hawkes_0cp.bug"
  )

  print(paste("Listo modelo", i, "-", d_names[i]))
  save(samples_hawkes0, file = "hawkes_0cp.RData")
}

Sys.time() - t1

# -----------------------------------------------------------------------------
# Quick inspection — change mod to inspect any of the 10 datasets
# -----------------------------------------------------------------------------
load("hawkes_0cp.RData")

mod <- 1
d       <- di[[mod]]
samples <- samples_hawkes0[[mod]]$BUGSoutput
y_name  <- d_names[mod]

su <- samples$summary

# Posterior summaries
su[grep("^alpha",  row.names(su)), ]
su[grep("^beta",   row.names(su)), ]
su[grep("^n$",     row.names(su)), ]
su[grep("^delta$", row.names(su)), ]
su[grep("^eta$",   row.names(su)), ]

# Branching ratio n: if posterior concentrates away from 0,
# autoexcitation is contributing meaningfully
cat("\nBranching ratio n (mean, sd):\n")
cat(su["n", c("mean", "sd")], "\n")

# -----------------------------------------------------------------------------
# Cumulative mean M(t) — observed vs estimated
# Note: in Hawkes, M(t) has no closed form because lambda depends on history.
# We use the Weibull baseline component as a reference, and compare
# cumulative observed exceedances with the posterior mean of A[k] to
# assess the clustering contribution.
# -----------------------------------------------------------------------------
acum <- sapply(1:t, function(i) sum(d <= i))

# Weibull baseline cumulative mean (same as NHPP model)
alpha_mean <- su["alpha", "mean"]
beta_mean  <- su["beta",  "mean"]
m_baseline <- (1:t / beta_mean)^alpha_mean

plot(1:t, acum,
     type = "l", lty = "dashed", lwd = 2,
     xlab = "Days", ylab = y_name,
     main = paste("Hawkes M1-0CP:", y_name))
lines(1:t, m_baseline, lwd = 2)
legend("topleft",
       legend = c("Observed cumulative", "Weibull baseline M(t)"),
       lty    = c("dashed", "solid"),
       lwd    = 2)

# -----------------------------------------------------------------------------
# Convergence diagnostics
# -----------------------------------------------------------------------------
source("07_utilidades/crear_mcmc_list.R")

# alpha
alphas <- samples$sims.list$alpha
alphas <- crear_mcmc_list(as.matrix(alphas), samples$n.chains)
densityplot(alphas, ylab = "", main = "Density of alpha")
xyplot(alphas, main = "Trace of alpha")
autocorr.plot(alphas)
gelman.diag(alphas)[1]
geweke.diag(alphas)
raftery.diag(alphas)
heidel.diag(alphas)

# beta
betas <- samples$sims.list$beta
betas <- crear_mcmc_list(as.matrix(betas), samples$n.chains)
densityplot(betas, ylab = "", main = "Density of beta")
xyplot(betas, main = "Trace of beta")
autocorr.plot(betas)
gelman.diag(betas)[1]
geweke.diag(betas)

# n (branching ratio)
ns <- samples$sims.list$n
ns <- crear_mcmc_list(as.matrix(ns), samples$n.chains)
densityplot(ns, ylab = "", main = "Density of n (branching ratio)")
xyplot(ns, main = "Trace of n")
autocorr.plot(ns)
gelman.diag(ns)[1]
geweke.diag(ns)

# delta
deltas <- samples$sims.list$delta
deltas <- crear_mcmc_list(as.matrix(deltas), samples$n.chains)
densityplot(deltas, ylab = "", main = "Density of delta")
xyplot(deltas, main = "Trace of delta")
autocorr.plot(deltas)
gelman.diag(deltas)[1]
geweke.diag(deltas)

# -----------------------------------------------------------------------------
# Model selection criteria
# -----------------------------------------------------------------------------

# DIC
samples$DIC

# WAIC — uses phi2 (log-likelihood contributions per observation)
loglik    <- samples$sims.list$phi2
waic_hawkes0 <- waic(loglik)
waic_hawkes0
loo(loglik)

# SDM (sum of absolute deviations between observed and estimated cumulative mean)
# Compare with NHPP SDM from published article
# sum(abs(acum[d] - m_baseline[d]))
