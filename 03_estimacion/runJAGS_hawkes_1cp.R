# =============================================================================
# Hawkes M1 - 1 Change Point
# Baseline: Weibull por segmentos (2 segmentos)
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
# Article values for 1CP: 80,000 iter, burnin 20,000, thin 20
# -----------------------------------------------------------------------------
niter   <- 80000
nburnin <- 20000
nthin   <- 20

# Constant C for zeros trick
C <- 10

# -----------------------------------------------------------------------------
# Estimation loop
# Uncomment the line below before first run:
samples_hawkes1 <- list()
# To resume after interruption, comment the line above, load the .RData,
# and change 1:4 to the remaining indices, e.g. 3:4
# load("hawkes_1cp.RData")
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
  
  samples_hawkes1[[i]] <- jags(
    data               = line_data,
    parameters.to.save = c("alpha", "beta", "n", "delta", "eta", "tau", "phi2"),
    n.chains           = 2,
    n.iter             = niter,
    n.burnin           = nburnin,
    n.thin             = nthin,
    model.file         = "02_modelos_JAGS/bug_files/hawkes_1cp.bug"
  )
  
  print(paste("Listo modelo", i, "-", d_names[i]))
  save(samples_hawkes1, file = "hawkes_1cp.RData")
}

Sys.time() - t1

# -----------------------------------------------------------------------------
# Inspection — runs over all 4 datasets
# -----------------------------------------------------------------------------
load("hawkes_1cp.RData")

source("07_utilidades/crear_mcmc_list.R")

for (mod in 1:4) {
  
  d       <- di[[mod]]
  samples <- samples_hawkes1[[mod]]$BUGSoutput
  y_name  <- d_names[mod]
  su      <- samples$summary
  
  cat(sprintf("\n\n====== Modelo %d | %s ======\n", mod, y_name))
  
  # Posterior summaries
  cat("\nalpha:\n"); print(su[grep("^alpha",  row.names(su)), c("mean","sd","2.5%","97.5%")])
  cat("\nbeta:\n");  print(su[grep("^beta",   row.names(su)), c("mean","sd","2.5%","97.5%")])
  cat("\ntau:\n");   print(su[grep("^tau",    row.names(su)), c("mean","sd","2.5%","97.5%")])
  cat("\nn:\n");     print(su[grep("^n$",     row.names(su)), c("mean","sd","2.5%","97.5%")])
  cat("\ndelta:\n"); print(su[grep("^delta$", row.names(su)), c("mean","sd","2.5%","97.5%")])
  cat("\neta:\n");   print(su[grep("^eta$",   row.names(su)), c("mean","sd","2.5%","97.5%")])
  
  cat(sprintf("\nBranching ratio n  — mean: %.4f  sd: %.4f\n",
              su["n",       "mean"], su["n",       "sd"]))
  cat(sprintf("Kernel decay delta — mean: %.4f  sd: %.4f\n",
              su["delta",   "mean"], su["delta",   "sd"]))
  cat(sprintf("Vida media kernel  — %.2f dias\n",
              log(2) / su["delta", "mean"]))
  cat(sprintf("tau[1]             — mean: %.1f  sd: %.2f\n",
              su[grep("^tau",   row.names(su)), "mean"], su[grep("^tau",   row.names(su)), "sd"]))
  
  # ---------------------------------------------------------------------------
  # Cumulative mean M(t) — Weibull baseline por segmentos
  # Segment 1: (t/beta[1])^alpha[1]                       t < tau[1]
  # Segment 2: m(tau[1]) + (t/beta[2])^alpha[2]
  #                      - (tau[1]/beta[2])^alpha[2]       t >= tau[1]
  # ---------------------------------------------------------------------------
  acum <- sapply(1:t, function(i) sum(d <= i))
  
  alpha_mu <- su[grep("^alpha", row.names(su)), "mean"]
  beta_mu  <- su[grep("^beta",  row.names(su)), "mean"]
  tau_mu   <- su[grep("^tau",   row.names(su)), "mean"]
  n_mu     <- su[grep("^n$",    row.names(su)), "mean"]
  delta_mu <- su[grep("^delta$",row.names(su)), "mean"]
  
  alpha1 <- alpha_mu[1];  beta1 <- beta_mu[1]
  alpha2 <- alpha_mu[2];  beta2 <- beta_mu[2]
  tau1   <- tau_mu[1]
  
  
  m_baseline <- ifelse(
    1:t < tau1,
    (1:t / beta1)^alpha1,
    (tau1 / beta1)^alpha1 + (1:t / beta2)^alpha2 - (tau1 / beta2)^alpha2
  )
  
  plot(1:t, acum,
       type = "l", lty = "dashed", lwd = 2,
       xlab = "Days", ylab = y_name,
       main = paste("Hawkes M1-1CP:", y_name))
  lines(1:t, m_baseline, lwd = 2)
  abline(v = tau1, col = "grey50", lty = 3)
  legend("topleft",
         legend = c("Observed cumulative", "Weibull baseline M(t)", "tau[1]"),
         lty    = c("dashed", "solid", "dotted"),
         lwd    = 2)
  
  cat(sprintf("\nSDM: %.4f\n", sum(abs(acum[d] - m_baseline[d]))))
  
  # ---------------------------------------------------------------------------
  # Convergence diagnostics
  # ---------------------------------------------------------------------------
  alphas <- crear_mcmc_list(as.matrix(samples$sims.list$alpha), samples$n.chains)
  betas  <- crear_mcmc_list(as.matrix(samples$sims.list$beta),  samples$n.chains)
  taus   <- crear_mcmc_list(as.matrix(samples$sims.list$tau),   samples$n.chains)
  ns     <- crear_mcmc_list(as.matrix(samples$sims.list$n),     samples$n.chains)
  deltas <- crear_mcmc_list(as.matrix(samples$sims.list$delta), samples$n.chains)
  
  # Alpha
  densityplot(alphas, ylab = "", main = paste("Density of alpha |", y_name))
  xyplot(alphas,      main = paste("Trace of alpha |",   y_name))
  autocorr.plot(alphas)
  print(gelman.diag(alphas)[1])
  print(geweke.diag(alphas))
  print(raftery.diag(alphas))
  print(heidel.diag(alphas))
  
  # Beta
  densityplot(betas, ylab = "", main = paste("Density of beta |", y_name))
  xyplot(betas,      main = paste("Trace of beta |",  y_name))
  autocorr.plot(betas)
  print(gelman.diag(betas)[1])
  print(geweke.diag(betas))
  
  # Tau
  densityplot(taus, ylab = "", main = paste("Density of tau |", y_name))
  xyplot(taus,      main = paste("Trace of tau |",  y_name))
  autocorr.plot(taus)
  print(gelman.diag(taus)[1])
  print(geweke.diag(taus))
  
  # n (branching ratio)
  densityplot(ns, ylab = "", main = paste("Density of n |", y_name))
  xyplot(ns,      main = paste("Trace of n |",  y_name))
  autocorr.plot(ns)
  print(gelman.diag(ns)[1])
  print(geweke.diag(ns))
  
  # delta
  densityplot(deltas, ylab = "", main = paste("Density of delta |", y_name))
  xyplot(deltas,      main = paste("Trace of delta |",  y_name))
  autocorr.plot(deltas)
  print(gelman.diag(deltas)[1])
  print(geweke.diag(deltas))
  
  # ---------------------------------------------------------------------------
  # Model selection criteria
  # ---------------------------------------------------------------------------
  cat("\nDIC:", samples$DIC, "\n")
  
  loglik <- samples$sims.list$phi2
  cat("\nWAIC:\n")
  print(waic(loglik))
  print(loo(loglik))
}