# Hawkes M1 - 2 Change Points

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
T_total <- 1096   # Días totales de la serie completa (referencia, no modificar)
T_train <- 916    # AJUSTAR AQUÍ: días usados para ajustar el modelo
#     Los días (T_train+1):T_total se reservan para predicción

di <- list(excePM25col_Bog, excePM25col_Med,
           excePM25who_Bog, excePM25who_Med)

d_names <- c("PM25col_Bog", "PM25col_Med",
             "PM25who_Bog", "PM25who_Med")

# -----------------------------------------------------------------------------
# MCMC settings
# -----------------------------------------------------------------------------
niter   <- 80000
nburnin <- 20000
nthin   <- 20

C <- 10

# -----------------------------------------------------------------------------
# Estimation loop
# -----------------------------------------------------------------------------
samples_hawkes2_forecast <- list()

t1 <- Sys.time()

for (i in 1:4) {
  
  d_full <- di[[i]]                      # serie completa de excedencias
  d      <- d_full[d_full <= T_train]    # solo el período de ajuste → JAGS
  K      <- length(d)
  
  line_data <- list(
    "K" = K,
    "T" = T_train,   # horizonte de integración = fin del período de ajuste
    "d" = d,
    "C" = C
  )
  
  samples_hawkes2_forecast[[i]] <- jags(
    data               = line_data,
    parameters.to.save = c("alpha", "beta", "n", "delta", "eta", "tau", "phi2"),
    n.chains           = 2,
    n.iter             = niter,
    n.burnin           = nburnin,
    n.thin             = nthin,
    model.file         = "02_modelos_JAGS/bug_files/hawkes_2cp.bug"
  )
  
  print(paste("Listo modelo", i, "-", d_names[i]))
  save(samples_hawkes2_forecast, file = "hawkes_2cp_train.RData")
}

Sys.time() - t1

# -----------------------------------------------------------------------------
# Quick inspection — cambiar mod para inspeccionar cada dataset
# -----------------------------------------------------------------------------
load("hawkes_2cp_train.RData")

mod <- 1
d_full  <- di[[mod]]
d       <- d_full[d_full <= T_train]
samples <- samples_hawkes2_forecast[[mod]]$BUGSoutput
y_name  <- d_names[mod]

su <- samples$summary

# Posterior summaries
su[grep("^alpha",  row.names(su)), ]
su[grep("^beta",   row.names(su)), ]
su[grep("^tau",    row.names(su)), ]
su[grep("^n$",     row.names(su)), ]
su[grep("^delta$", row.names(su)), ]
su[grep("^eta$",   row.names(su)), ]

cat("\nBranching ratio n (mean, sd):\n")
cat(su["n", c("mean", "sd")], "\n")

cat("\nChange points tau (mean, sd):\n")
print(su[grep("^tau", row.names(su)), c("mean", "sd")])

# -----------------------------------------------------------------------------
# M(t): observado vs baseline Weibull — período de AJUSTE
# Segment 1: (t/beta[1])^alpha[1]                               t < tau[1]
# Segment 2: m(tau[1]) + (t/beta[2])^alpha[2]
#                      - (tau[1]/beta[2])^alpha[2]    tau[1] <= t < tau[2]
# Segment 3: m(tau[2]) + (t/beta[3])^alpha[3]
#                      - (tau[2]/beta[3])^alpha[3]               t >= tau[2]
# -----------------------------------------------------------------------------
acum_train <- sapply(1:T_train, function(i) sum(d <= i))

alpha1 <- su["alpha[1]", "mean"];  beta1 <- su["beta[1]", "mean"]
alpha2 <- su["alpha[2]", "mean"];  beta2 <- su["beta[2]", "mean"]
alpha3 <- su["alpha[3]", "mean"];  beta3 <- su["beta[3]", "mean"]
tau1   <- su["tau[1]",   "mean"]
tau2   <- su["tau[2]",   "mean"]

m_tau1 <- (tau1 / beta1)^alpha1
m_tau2 <- m_tau1 + (tau2 / beta2)^alpha2 - (tau1 / beta2)^alpha2

m_baseline <- ifelse(
  1:T_train < tau1,
  (1:T_train / beta1)^alpha1,
  ifelse(
    1:T_train < tau2,
    m_tau1 + (1:T_train / beta2)^alpha2 - (tau1 / beta2)^alpha2,
    m_tau2 + (1:T_train / beta3)^alpha3 - (tau2 / beta3)^alpha3
  )
)

plot(1:T_train, acum_train,
     type = "l", lty = "dashed", lwd = 2,
     xlab = "Days", ylab = y_name,
     main = paste("Hawkes M1-2CP (ajuste, T_train =", T_train, "):", y_name))
lines(1:T_train, m_baseline, lwd = 2)
abline(v = c(tau1, tau2), col = "grey50", lty = 3)
legend("topleft",
       legend = c("Observed cumulative (train)", "Weibull baseline M(t)", "tau[1], tau[2]"),
       lty    = c("dashed", "solid", "dotted"),
       lwd    = 2)

# -----------------------------------------------------------------------------
# Comparación predictiva: período (T_train+1):T_total
# -----------------------------------------------------------------------------
d_test     <- d_full[d_full > T_train]
n_test_obs <- length(d_test)

# Conteo acumulado sobre la serie completa
acum_full <- sapply(1:T_total, function(i) sum(d_full <= i))

# Extrapolación del baseline Weibull al período de predicción
m_baseline_full <- ifelse(
  1:T_total < tau1,
  (1:T_total / beta1)^alpha1,
  ifelse(
    1:T_total < tau2,
    m_tau1 + (1:T_total / beta2)^alpha2 - (tau1 / beta2)^alpha2,
    m_tau2 + (1:T_total / beta3)^alpha3 - (tau2 / beta3)^alpha3
  )
)

pred_expected <- m_baseline_full[T_total] - m_baseline_full[T_train]

cat("\n--- Comparación predictiva ---\n")
cat("Período de predicción: días", T_train + 1, "a", T_total, "\n")
cat("Excedencias observadas en hold-out:", n_test_obs, "\n")
cat("Excedencias esperadas (baseline Weibull):", round(pred_expected, 2), "\n")

# Gráfica: serie completa con línea vertical en T_train
plot(1:T_total, acum_full,
     type = "l", lty = "dashed", lwd = 2, col = "black",
     xlab = "Days", ylab = y_name,
     main = paste("Hawkes M1-2CP (ajuste + predicción):", y_name))
lines(1:T_total, m_baseline_full, lwd = 2, col = "blue")
abline(v = c(tau1, tau2), col = "grey50", lty = 3)
abline(v = T_train, col = "red", lty = 2, lwd = 1.5)
text(T_train + 10, max(acum_full) * 0.5,
     paste("T_train =", T_train), col = "red", cex = 0.8, adj = 0)
legend("topleft",
       legend = c("Observed cumulative (full)",
                  "Weibull baseline M(t) (extrapolated)",
                  "tau[1], tau[2]",
                  "End of training period"),
       col    = c("black", "blue", "grey50", "red"),
       lty    = c("dashed", "solid", "dotted", "dashed"),
       lwd    = 2)

# -----------------------------------------------------------------------------
# Convergence diagnostics
# -----------------------------------------------------------------------------
source("07_utilidades/crear_mcmc_list.R")

alphas <- samples$sims.list$alpha
alphas <- crear_mcmc_list(as.matrix(alphas), samples$n.chains)
densityplot(alphas, ylab = "", main = "Density of alpha")
xyplot(alphas, main = "Trace of alpha")
autocorr.plot(alphas)
gelman.diag(alphas)[1]
geweke.diag(alphas)
raftery.diag(alphas)
heidel.diag(alphas)

betas <- samples$sims.list$beta
betas <- crear_mcmc_list(as.matrix(betas), samples$n.chains)
densityplot(betas, ylab = "", main = "Density of beta")
xyplot(betas, main = "Trace of beta")
autocorr.plot(betas)
gelman.diag(betas)[1]
geweke.diag(betas)

taus <- samples$sims.list$tau
taus <- crear_mcmc_list(as.matrix(taus), samples$n.chains)
densityplot(taus, ylab = "", main = "Density of tau")
xyplot(taus, main = "Trace of tau")
autocorr.plot(taus)
gelman.diag(taus)[1]
geweke.diag(taus)

ns <- samples$sims.list$n
ns <- crear_mcmc_list(as.matrix(ns), samples$n.chains)
densityplot(ns, ylab = "", main = "Density of n (branching ratio)")
xyplot(ns, main = "Trace of n")
autocorr.plot(ns)
gelman.diag(ns)[1]
geweke.diag(ns)

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
samples$DIC

loglik              <- samples$sims.list$phi2
waic_hawkes2_train  <- waic(loglik)
waic_hawkes2_train
loo(loglik)