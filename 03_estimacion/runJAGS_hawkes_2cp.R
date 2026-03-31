# =============================================================================
# Hawkes M1 - 2 Change Points
# Baseline: Weibull por segmentos (3 segmentos)
# Kernel:   exponential, eta and delta global
# Model:    lambda(t) = (alpha_j/beta_j)*(t/beta_j)^(alpha_j-1) + eta * A(t)
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
# 2CP: same as NHPP-2CP from article (80000 iter, burnin 20000, thin 20)
# -----------------------------------------------------------------------------
niter   <- 80000
nburnin <- 20000
nthin   <- 20

# Constant C for zeros trick: must satisfy phi[k] = -log(lambda_k) + l2 + l3 + C >= 0
# C = 10 is conservative for typical exceedance intensities
C <- 10

# -----------------------------------------------------------------------------
# Estimation loop
# Uncomment the line below before first run:
samples_hawkes2 <- list()
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

  samples_hawkes2[[i]] <- jags(
    data               = line_data,
    parameters.to.save = c("alpha", "beta", "n", "delta", "eta", "tau", "phi2"),
    n.chains           = 2,
    n.iter             = niter,
    n.burnin           = nburnin,
    n.thin             = nthin,
    model.file         = "02_modelos_JAGS/bug_files/hawkes_2cp.bug"
  )

  print(paste("Listo modelo", i, "-", d_names[i]))
  save(samples_hawkes2, file = "hawkes_2cp.RData")
}

Sys.time() - t1

# -----------------------------------------------------------------------------
# Quick inspection — change mod to inspect any of the 4 datasets
# -----------------------------------------------------------------------------
load("hawkes_2cp.RData")

mod <- 1
d       <- di[[mod]]
samples <- samples_hawkes2[[mod]]$BUGSoutput
y_name  <- d_names[mod]

su <- samples$summary

# Posterior summaries
su[grep("^alpha",  row.names(su)), ]
su[grep("^beta",   row.names(su)), ]
su[grep("^tau",    row.names(su)), ]
su[grep("^n$",     row.names(su)), ]
su[grep("^delta$", row.names(su)), ]
su[grep("^eta$",   row.names(su)), ]

# Branching ratio n: if posterior concentrates away from 0,
# autoexcitation is contributing meaningfully
cat("\nBranching ratio n (mean, sd):\n")
cat(su["n", c("mean", "sd")], "\n")

cat("\nChange points tau (mean, sd):\n")
print(su[grep("^tau", row.names(su)), c("mean", "sd")])

# -----------------------------------------------------------------------------
# Cumulative mean M(t) — observed vs estimated
# Segment 1: (t/beta[1])^alpha[1]                               t < tau[1]
# Segment 2: m(tau[1]) + (t/beta[2])^alpha[2]
#                      - (tau[1]/beta[2])^alpha[2]    tau[1] <= t < tau[2]
# Segment 3: m(tau[2]) + (t/beta[3])^alpha[3]
#                      - (tau[2]/beta[3])^alpha[3]               t >= tau[2]
# -----------------------------------------------------------------------------
acum <- sapply(1:t, function(i) sum(d <= i))

alpha1 <- su["alpha[1]", "mean"];  beta1 <- su["beta[1]", "mean"]
alpha2 <- su["alpha[2]", "mean"];  beta2 <- su["beta[2]", "mean"]
alpha3 <- su["alpha[3]", "mean"];  beta3 <- su["beta[3]", "mean"]
tau1   <- su["tau[1]",   "mean"]
tau2   <- su["tau[2]",   "mean"]

m_tau1 <- (tau1 / beta1)^alpha1
m_tau2 <- m_tau1 + (tau2 / beta2)^alpha2 - (tau1 / beta2)^alpha2

m_baseline <- ifelse(
  1:t < tau1,
  (1:t / beta1)^alpha1,
  ifelse(
    1:t < tau2,
    m_tau1 + (1:t / beta2)^alpha2 - (tau1 / beta2)^alpha2,
    m_tau2 + (1:t / beta3)^alpha3 - (tau2 / beta3)^alpha3
  )
)

plot(1:t, acum,
     type = "l", lty = "dashed", lwd = 2,
     xlab = "Days", ylab = y_name,
     main = paste("Hawkes M1-2CP:", y_name))
lines(1:t, m_baseline, lwd = 2)
abline(v = c(tau1, tau2), col = "grey50", lty = 3)
legend("topleft",
       legend = c("Observed cumulative", "Weibull baseline M(t)", "tau[1], tau[2]"),
       lty    = c("dashed", "solid", "dotted"),
       lwd    = 2)

# -----------------------------------------------------------------------------
# Convergence diagnostics
# -----------------------------------------------------------------------------
source("07_utilidades/crear_mcmc_list.R")

# alpha (matrix: iterations x 3 segments)
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

# tau
taus <- samples$sims.list$tau
taus <- crear_mcmc_list(as.matrix(taus), samples$n.chains)
densityplot(taus, ylab = "", main = "Density of tau")
xyplot(taus, main = "Trace of tau")
autocorr.plot(taus)
gelman.diag(taus)[1]
geweke.diag(taus)

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
loglik       <- samples$sims.list$phi2
waic_hawkes2 <- waic(loglik)
waic_hawkes2
loo(loglik)

# SDM (sum of absolute deviations between observed and estimated cumulative mean)
# Compare with NHPP SDM from published article
# sum(abs(acum[d] - m_baseline[d]))
