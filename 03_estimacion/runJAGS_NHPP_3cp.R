# =============================================================================
# NHPP Weibull - 3 Change Points
# Baseline: Weibull por segmentos (4 segmentos)
# tau[1]~dunif(0,365), tau[2]~dunif(365,730), tau[3]~dunif(730,1096)
# beta[i] ~ dunif(0,100) para PM2.5
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
t <- 1096

di <- list(excePM25col_Bog, excePM25col_Med,
           excePM25who_Bog, excePM25who_Med)

d_names <- c("PM25col_Bog", "PM25col_Med",
             "PM25who_Bog", "PM25who_Med")

# -----------------------------------------------------------------------------
# MCMC settings — articulo: 100,000 iter, burnin 20,000, thin 20
# -----------------------------------------------------------------------------
niter   <- 100000
nburnin <- 20000
nthin   <- 20

# inits: una lista completa por cadena con alpha y beta juntos
line_inits <- list(
  list(alpha = c(1, 1, 1, 1), beta = c(1, 1, 1, 1)),
  list(alpha = c(1, 1, 1, 1), beta = c(1, 1, 1, 1))
)

# -----------------------------------------------------------------------------
# Estimation loop
# Uncomment before first run:
samples_NHPP3 <- list()
# To resume: comment line above, load RData, change 1:4 to remaining indices
# load("NHPP_3cp.RData")
# -----------------------------------------------------------------------------

t1 <- Sys.time()

for (i in 1:4) {
  
  d <- di[[i]]
  K <- length(d)
  
  line_data <- list("K" = K, "T" = t, "d" = d)
  
  samples_NHPP3[[i]] <- jags(
    data               = line_data,
    inits              = line_inits,
    parameters.to.save = c("alpha", "beta", "m", "tau", "phi2"),
    n.chains           = 2,
    n.iter             = niter,
    n.burnin           = nburnin,
    n.thin             = nthin,
    model.file         = "02_modelos_JAGS/bug_files/NHPP_3cp.bug"
  )
  
  print(paste("Listo modelo", i, "-", d_names[i]))
  save(samples_NHPP3, file = "NHPP_3cp.RData")
}

Sys.time() - t1

# -----------------------------------------------------------------------------
# Inspection — all 4 datasets
# -----------------------------------------------------------------------------
load("NHPP_3cp.RData")

source("07_utilidades/crear_mcmc_list.R")

for (mod in 1:4) {
  
  d       <- di[[mod]]
  samples <- samples_NHPP3[[mod]]$BUGSoutput
  y_name  <- d_names[mod]
  su      <- samples$summary
  
  cat(sprintf("\n\n====== Modelo %d | %s ======\n", mod, y_name))
  
  # Posterior summaries
  cat("\nalpha:\n"); print(su[grep("^alpha", row.names(su)), c("mean","sd","2.5%","97.5%")])
  cat("\nbeta:\n");  print(su[grep("^beta",  row.names(su)), c("mean","sd","2.5%","97.5%")])
  cat("\ntau:\n");   print(su[grep("^tau",   row.names(su)), c("mean","sd","2.5%","97.5%")])
  
  tau_su <- su[grep("^tau", row.names(su)), , drop = FALSE]
  cat(sprintf("\ntau[1] — mean: %.1f  sd: %.2f\n", tau_su[1,"mean"], tau_su[1,"sd"]))
  cat(sprintf("tau[2] — mean: %.1f  sd: %.2f\n",   tau_su[2,"mean"], tau_su[2,"sd"]))
  cat(sprintf("tau[3] — mean: %.1f  sd: %.2f\n",   tau_su[3,"mean"], tau_su[3,"sd"]))
  
  # ---------------------------------------------------------------------------
  # Cumulative mean M(t) — observed vs estimated
  # m[i] evaluated at exceedance days d[i] — from .bug
  # ---------------------------------------------------------------------------
  acum  <- sapply(1:t, function(i) sum(d <= i))
  mm    <- su[grep("^m", row.names(su)), "mean"]
  mm_lo <- su[grep("^m", row.names(su)), "2.5%"]
  mm_hi <- su[grep("^m", row.names(su)), "97.5%"]
  
  tau1 <- tau_su[1, "mean"]
  tau2 <- tau_su[2, "mean"]
  tau3 <- tau_su[3, "mean"]
  
  plot(1:t, acum,
       type = "l", lty = "dashed", lwd = 2,
       xlab = "Days", ylab = y_name,
       main = paste("NHPP-3CP:", y_name),
       ylim = c(0, max(c(acum, mm_hi), na.rm = TRUE)))
  lines(d, mm,    col = "blue", lwd = 2)
  lines(d, mm_lo, col = "red",  lty = "dashed")
  lines(d, mm_hi, col = "red",  lty = "dashed")
  abline(v = c(tau1, tau2, tau3), col = "grey50", lty = 3)
  legend("topleft",
         legend = c("Observed cumulative", "NHPP estimated M(t)", "95% CI", "tau[1,2,3]"),
         lty    = c("dashed", "solid", "dashed", "dotted"),
         col    = c("black", "blue", "red", "grey50"), lwd = 2)
  
  cat(sprintf("\nSDM: %.4f\n", sum(abs(acum[d] - mm))))
  
  # ---------------------------------------------------------------------------
  # Lambda por segmentos — plot
  # ---------------------------------------------------------------------------
  alpha_mu <- su[grep("^alpha", row.names(su)), "mean"]
  beta_mu  <- su[grep("^beta",  row.names(su)), "mean"]
  
  lambda <- ifelse(d < tau1,
                   (alpha_mu[1]/beta_mu[1]) * (d/beta_mu[1])^(alpha_mu[1]-1),
                   ifelse(d < tau2,
                          (alpha_mu[2]/beta_mu[2]) * (d/beta_mu[2])^(alpha_mu[2]-1),
                          ifelse(d < tau3,
                                 (alpha_mu[3]/beta_mu[3]) * (d/beta_mu[3])^(alpha_mu[3]-1),
                                 (alpha_mu[4]/beta_mu[4]) * (d/beta_mu[4])^(alpha_mu[4]-1)
                          )
                   )
  )
  
  sep <- sapply(c(tau1, tau2, tau3), function(tk) which(d >= tk)[1])
  sep <- c(1, sep, length(d) + 1)
  
  plot(d, lambda, type = "n", xlab = "Days", ylab = y_name,
       main = paste("Lambda NHPP-3CP:", y_name))
  for (s in 1:(length(sep)-1)) {
    idx <- sep[s]:(sep[s+1]-1)
    lines(d[idx], lambda[idx], lty = "solid", lwd = 1.5)
  }
  abline(v = c(tau1, tau2, tau3), col = "grey50", lty = 3)
  
  # ---------------------------------------------------------------------------
  # Convergence diagnostics
  # ---------------------------------------------------------------------------
  alphas <- crear_mcmc_list(as.matrix(samples$sims.list$alpha), samples$n.chains)
  betas  <- crear_mcmc_list(as.matrix(samples$sims.list$beta),  samples$n.chains)
  taus   <- crear_mcmc_list(as.matrix(samples$sims.list$tau),   samples$n.chains)
  
  # Alpha
  densityplot(alphas, ylab = "", main = paste("Density alpha |", y_name))
  xyplot(alphas,      main = paste("Trace alpha |",   y_name))
  autocorr.plot(alphas)
  print(gelman.diag(alphas)[1])
  print(geweke.diag(alphas))
  print(raftery.diag(alphas))
  print(heidel.diag(alphas))
  
  # Beta
  densityplot(betas, ylab = "", main = paste("Density beta |", y_name))
  xyplot(betas,      main = paste("Trace beta |",  y_name))
  autocorr.plot(betas)
  print(gelman.diag(betas)[1])
  print(geweke.diag(betas))
  
  # Tau
  densityplot(taus, ylab = "", main = paste("Density tau |", y_name))
  xyplot(taus,      main = paste("Trace tau |",  y_name))
  autocorr.plot(taus)
  print(gelman.diag(taus)[1])
  print(geweke.diag(taus))
  
  # ---------------------------------------------------------------------------
  # Model selection criteria
  # ---------------------------------------------------------------------------
  cat("\nDIC:", round(samples$DIC, 0), "\n")
  
  loglik <- samples$sims.list$phi2
  cat("\nWAIC:\n")
  print(waic(loglik))
  print(loo(loglik))
}