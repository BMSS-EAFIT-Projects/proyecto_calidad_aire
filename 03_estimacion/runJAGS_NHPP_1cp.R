# =============================================================================
# NHPP Weibull - 1 Change Point
# Baseline: Weibull por segmentos (2 segmentos)
# Zeros trick: phi[i] = -log(lambda[i] * F),  F absorbs integral over [0,T]
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

# -----------------------------------------------------------------------------
# Estimation loop
# Uncomment the line below before first run:
samples_NHPP1 <- list()
# -----------------------------------------------------------------------------

t1 <- Sys.time()

for (i in 1:4) {
  
  d <- di[[i]]
  K <- length(d)
  
  line_data <- list(
    "K" = K,
    "T" = t,
    "d" = d
  )
  
  samples_NHPP1[[i]] <- jags(
    data               = line_data,
    parameters.to.save = c("alpha", "beta", "m", "tau", "phi2"),
    n.chains           = 2,
    n.iter             = niter,
    n.burnin           = nburnin,
    n.thin             = nthin,
    model.file         = "02_modelos_JAGS/bug_files/NHPP_1cp.bug"
  )
  
  print(paste("Listo modelo", i, "-", d_names[i]))
  save(samples_NHPP1, file = "NHPP_1cp.RData")
}

Sys.time() - t1

# -----------------------------------------------------------------------------
# Inspección — recorre los 4 contaminantes
# -----------------------------------------------------------------------------
load("NHPP_1cp.RData")

source("07_utilidades/crear_mcmc_list.R")

for (mod in 1:4) {
  
  d       <- di[[mod]]
  samples <- samples_NHPP1[[mod]]$BUGSoutput
  y_name  <- d_names[mod]
  su      <- samples$summary
  
  cat("\n", rep("=", 60), "\n")
  cat("Modelo:", mod, "|", y_name, "\n")
  cat(rep("=", 60), "\n")
  
  # Posterior summaries
  cat("\nalpha:\n");  print(su[grep("^alpha", row.names(su)), c("mean","sd","2.5%","97.5%")])
  cat("\nbeta:\n");   print(su[grep("^beta",  row.names(su)), c("mean","sd","2.5%","97.5%")])
  cat("\ntau:\n");    print(su[grep("^tau",   row.names(su)), c("mean","sd","2.5%","97.5%")])
  
  # ---------------------------------------------------------------------------
  # Cumulative mean M(t)
  # ---------------------------------------------------------------------------
  acum  <- sapply(1:t, function(i) sum(d <= i))
  mm    <- su[grep("^m", row.names(su)), "mean"]
  mm_lo <- su[grep("^m", row.names(su)), "2.5%"]
  mm_hi <- su[grep("^m", row.names(su)), "97.5%"]
  tau1  <- tau1 <- su[grep("^tau", row.names(su)), "mean"][1]
  
  plot(1:t, acum,
       type = "l", lty = "dashed", lwd = 2,
       xlab = "Days", ylab = y_name,
       main = paste("NHPP-1CP:", y_name))
  lines(d, mm,    col = "blue", lwd = 2)
  lines(d, mm_lo, col = "red",  lty = "dashed")
  lines(d, mm_hi, col = "red",  lty = "dashed")
  abline(v = tau1, col = "grey50", lty = 3)
  legend("topleft",
         legend = c("Observed", "NHPP M(t)", "95% CI", "tau[1]"),
         lty    = c("dashed", "solid", "dashed", "dotted"),
         col    = c("black", "blue", "red", "grey50"),
         lwd    = 2)
  
  cat("\nSDM:", sum(abs(acum[d] - mm)), "\n")
  
  # ---------------------------------------------------------------------------
  # Convergence diagnostics
  # ---------------------------------------------------------------------------
  alphas <- crear_mcmc_list(as.matrix(samples$sims.list$alpha), samples$n.chains)
  betas  <- crear_mcmc_list(as.matrix(samples$sims.list$beta),  samples$n.chains)
  taus   <- crear_mcmc_list(as.matrix(samples$sims.list$tau),   samples$n.chains)
  
  densityplot(alphas, ylab = "", main = paste(y_name, "- Density alpha"))
  xyplot(alphas,      main = paste(y_name, "- Trace alpha"))
  autocorr.plot(alphas)
  cat("\nGelman alpha:\n");  print(gelman.diag(alphas)[1])
  cat("\nGeweke alpha:\n");  print(geweke.diag(alphas))
  raftery.diag(alphas)
  heidel.diag(alphas)
  
  densityplot(betas, ylab = "", main = paste(y_name, "- Density beta"))
  xyplot(betas,      main = paste(y_name, "- Trace beta"))
  autocorr.plot(betas)
  cat("\nGelman beta:\n");   print(gelman.diag(betas)[1])
  cat("\nGeweke beta:\n");   print(geweke.diag(betas))
  
  densityplot(taus, ylab = "", main = paste(y_name, "- Density tau"))
  xyplot(taus,      main = paste(y_name, "- Trace tau"))
  autocorr.plot(taus)
  cat("\nGelman tau:\n");    print(gelman.diag(taus)[1])
  cat("\nGeweke tau:\n");    print(geweke.diag(taus))
  
  # ---------------------------------------------------------------------------
  # Model selection criteria
  # ---------------------------------------------------------------------------
  cat("\nDIC:", samples$DIC, "\n")
  
  loglik <- samples$sims.list$phi2
  cat("\nWAIC:\n")
  print(waic(loglik))
}