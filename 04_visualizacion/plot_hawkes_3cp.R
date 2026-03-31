# =============================================================================
# Visualización Hawkes M4 - 3 Change Points  [CORREGIDO]
# =============================================================================
# CONVENCIÓN DE TIEMPO (corrección aplicada):
#   BUGS estima beta[k] con tiempo ABSOLUTO en todos los segmentos:
#     Seg 1: (tau[1]/beta[1])^alpha[1]
#     Seg 2: (tau[2]/beta[2])^alpha[2] - (tau[1]/beta[2])^alpha[2]
#     Seg 3: (tau[3]/beta[3])^alpha[3] - (tau[2]/beta[3])^alpha[3]
#     Seg 4: (T/beta[4])^alpha[4]      - (tau[3]/beta[4])^alpha[4]
#   En R se usa la misma convención:
#     M_k(t) = (t_end/beta[k])^alpha[k] - (lo[k]/beta[k])^alpha[k]
#   mu_k(t) = (alpha[k]/beta[k]) * (t/beta[k])^(alpha[k]-1)  [t absoluto]
# =============================================================================

.seg_lo_3cp <- function(tau) c(0, tau[1], tau[2], tau[3])
.seg_hi_3cp <- function(tau) c(tau[1], tau[2], tau[3], Inf)
.seg_k_3cp  <- function(ti, tau) {
  if      (ti < tau[1]) 1L
  else if (ti < tau[2]) 2L
  else if (ti < tau[3]) 3L
  else                  4L
}

# --------------------------------------------------------------------------
# Auxiliar: lambda(t | H_t) — 3 CPs
# --------------------------------------------------------------------------
calc_lambda_hawkes_3cp <- function(t_grid, d, tau,
                                   alpha, beta, n_par, delta) {
  eta <- n_par * delta
  sapply(t_grid, function(ti) {
    k <- .seg_k_3cp(ti, tau)
    if (ti <= 0) return(0)
    mu_t <- (alpha[k] / beta[k]) * (ti / beta[k])^(alpha[k] - 1)
    prev <- d[d < ti]
    kern <- if (length(prev) > 0) eta * sum(exp(-delta * (ti - prev))) else 0
    mu_t + kern
  })
}

# --------------------------------------------------------------------------
# Auxiliar: M(t) Hawkes exacta — 3 CPs  [CORREGIDO: tiempo absoluto]
# --------------------------------------------------------------------------
calc_M_hawkes_exacta_3cp <- function(t_eval, d, tau,
                                     alpha, beta, n_par, delta) {
  eta <- n_par * delta
  lo  <- .seg_lo_3cp(tau)
  hi  <- .seg_hi_3cp(tau)
  sapply(t_eval, function(tt) {
    M <- 0
    for (k in 1:4) {
      if (tt <= lo[k]) break
      t_end_k <- min(tt, hi[k])
      M <- M + (t_end_k / beta[k])^alpha[k] - (lo[k] / beta[k])^alpha[k]
      dk <- d[d >= lo[k] & d < t_end_k]
      if (length(dk) > 0)
        M <- M + (eta / delta) * sum(1 - exp(-delta * (tt - dk)))
    }
    M
  })
}

# --------------------------------------------------------------------------
# Auxiliar: bandas M(t) Hawkes — 3 CPs
# --------------------------------------------------------------------------
calc_M_hawkes_banda_3cp <- function(t_eval, d, tau_sims,
                                    alpha_sims, beta_sims,
                                    n_sims, delta_sims) {
  S        <- nrow(alpha_sims)
  m_matrix <- matrix(NA, nrow=S, ncol=length(t_eval))
  for (s in 1:S)
    m_matrix[s, ] <- calc_M_hawkes_exacta_3cp(
      t_eval, d, tau_sims[s, ],
      alpha_sims[s, ], beta_sims[s, ], n_sims[s], delta_sims[s])
  list(mean=apply(m_matrix,2,mean),
       lo  =apply(m_matrix,2,quantile,0.025),
       hi  =apply(m_matrix,2,quantile,0.975))
}

# --------------------------------------------------------------------------
# Auxiliar: M(t) NHPP exacta — 3 CPs  [CORREGIDO: tiempo absoluto]
# --------------------------------------------------------------------------
calc_M_nhpp_exacta_3cp <- function(t_eval, tau, alpha, beta) {
  lo <- .seg_lo_3cp(tau); hi <- .seg_hi_3cp(tau)
  sapply(t_eval, function(tt) {
    M <- 0
    for (k in 1:4) {
      if (tt <= lo[k]) break
      t_end_k <- min(tt, hi[k])
      M <- M + (t_end_k / beta[k])^alpha[k] - (lo[k] / beta[k])^alpha[k]
    }; M
  })
}

# --------------------------------------------------------------------------
# Auxiliar: bandas M(t) NHPP — 3 CPs  [CORREGIDO: tiempo absoluto]
# --------------------------------------------------------------------------
calc_M_nhpp_banda_3cp <- function(t_eval, tau_sims, alpha_sims, beta_sims) {
  S        <- nrow(alpha_sims)
  m_matrix <- matrix(NA, nrow=S, ncol=length(t_eval))
  for (s in 1:S) {
    a <- alpha_sims[s, ]; b <- beta_sims[s, ]; tau <- tau_sims[s, ]
    lo <- .seg_lo_3cp(tau); hi <- .seg_hi_3cp(tau)
    m_matrix[s, ] <- sapply(t_eval, function(tt) {
      M <- 0
      for (k in 1:4) {
        if (tt <= lo[k]) break
        t_end_k <- min(tt, hi[k])
        M <- M + (t_end_k / b[k])^a[k] - (lo[k] / b[k])^a[k]
      }; M
    })
  }
  list(mean=apply(m_matrix,2,mean),
       lo  =apply(m_matrix,2,quantile,0.025),
       hi  =apply(m_matrix,2,quantile,0.975))
}

# --------------------------------------------------------------------------
# Auxiliar: baseline Weibull — 3 CPs  [CORREGIDO: t absoluto]
# --------------------------------------------------------------------------
calc_mu_weibull_3cp <- function(t_grid, tau, alpha, beta) {
  sapply(t_grid, function(ti) {
    k <- .seg_k_3cp(ti, tau)
    if (ti <= 0) return(0)
    (alpha[k] / beta[k]) * (ti / beta[k])^(alpha[k] - 1)
  })
}

# =============================================================================
# 1. lambda(t | H_t) — 3 CPs
# =============================================================================
plot_lambda_hawkes_3cp <- function(n_modelo,
                                   samples_hawkes  = samples_hawkes3,
                                   di,
                                   y_name          = "Lambda",
                                   grilla_paso     = 5,
                                   mostrar_eventos = TRUE,
                                   col_lambda      = "black",
                                   col_mu          = "gray50",
                                   col_eventos     = "red",
                                   col_cp          = "purple",
                                   lwd_lambda      = 1.5) {
  d  <- di[[n_modelo]]
  su <- samples_hawkes[[n_modelo]]$BUGSoutput$summary
  alpha <- su[c("alpha[1]","alpha[2]","alpha[3]","alpha[4]"), "mean"]
  beta  <- su[c("beta[1]", "beta[2]", "beta[3]", "beta[4]"),  "mean"]
  n_par <- su["n",     "mean"]
  delta <- su["delta", "mean"]
  tau   <- su[c("tau[1]","tau[2]","tau[3]"), "mean"]
  cat(sprintf("\nModelo %d: %s — 3 CPs\n", n_modelo, y_name))
  cat(sprintf("  tau = [%.1f, %.1f, %.1f]\n", tau[1], tau[2], tau[3]))
  cat(sprintf("  n=%.4f  delta=%.4f  eta=%.4f  t1/2=%.2f\n",
              n_par, delta, n_par*delta, log(2)/delta))
  for (k in 1:4)
    cat(sprintf("  Seg%d: alpha=%.4f  beta=%.4f\n", k, alpha[k], beta[k]))
  t_grid      <- seq(1, 1096, by=grilla_paso)
  lambda_full <- calc_lambda_hawkes_3cp(t_grid, d, tau, alpha, beta, n_par, delta)
  mu_base     <- calc_mu_weibull_3cp(t_grid, tau, alpha, beta)
  plot(t_grid, lambda_full, type="l", lwd=lwd_lambda, col=col_lambda,
       xlab="Days", ylab=y_name, main=paste("Hawkes M4-3CP:", y_name),
       ylim=c(0, max(lambda_full, na.rm=TRUE)*1.1))
  lines(t_grid, mu_base, lwd=1, col=col_mu, lty="dashed")
  abline(v=tau, col=col_cp, lty="dotted", lwd=1.5)
  if (mostrar_eventos)
    rug(d, col=adjustcolor(col_eventos, alpha.f=0.4), ticksize=0.03)
  legend("topright", legend=c("lambda(t | H_t)","mu(t) baseline","CPs tau1/2/3"),
         col=c(col_lambda,col_mu,col_cp),
         lty=c("solid","dashed","dotted"), lwd=c(lwd_lambda,1,1.5), cex=0.8)
}

# =============================================================================
# 2. Baseline mu(t): NHPP vs Hawkes — 3 CPs
# =============================================================================
plot_baseline_comparison_3cp <- function(n_modelo,
                                         samples_nhpp   = samples_NHPP3,
                                         samples_hawkes = samples_hawkes3,
                                         y_name         = "mu(t)",
                                         col_nhpp       = "blue",
                                         col_hawkes     = "darkgreen",
                                         col_cp         = "purple") {
  t_grid  <- 1:1096
  su_nhpp <- samples_nhpp[[n_modelo]]$BUGSoutput$summary
  su_hwk  <- samples_hawkes[[n_modelo]]$BUGSoutput$summary
  alpha_nhpp <- su_nhpp[c("alpha[1]","alpha[2]","alpha[3]","alpha[4]"), "mean"]
  beta_nhpp  <- su_nhpp[c("beta[1]", "beta[2]", "beta[3]", "beta[4]"),  "mean"]
  tau_nhpp   <- su_nhpp[c("tau[1]","tau[2]","tau[3]"), "mean"]
  alpha_hwk  <- su_hwk[c("alpha[1]","alpha[2]","alpha[3]","alpha[4]"), "mean"]
  beta_hwk   <- su_hwk[c("beta[1]", "beta[2]", "beta[3]", "beta[4]"),  "mean"]
  tau_hwk    <- su_hwk[c("tau[1]","tau[2]","tau[3]"), "mean"]
  mu_nhpp <- calc_mu_weibull_3cp(t_grid, tau_nhpp, alpha_nhpp, beta_nhpp)
  mu_hwk  <- calc_mu_weibull_3cp(t_grid, tau_hwk,  alpha_hwk,  beta_hwk)
  plot(t_grid, mu_nhpp, type="l", lwd=2, col=col_nhpp,
       xlab="Days", ylab=y_name,
       main=paste("Baseline NHPP vs Hawkes 3CP:", y_name),
       ylim=c(0, max(c(mu_nhpp,mu_hwk), na.rm=TRUE)*1.1))
  lines(t_grid, mu_hwk, lwd=2, col=col_hawkes, lty="dashed")
  abline(v=tau_hwk, col=col_cp, lty="dotted", lwd=1.5)
  legend("topright", legend=c("NHPP","Hawkes M4","CPs tau1/2/3"),
         col=c(col_nhpp,col_hawkes,col_cp),
         lty=c("solid","dashed","dotted"), lwd=2, cex=0.8)
}

# =============================================================================
# 3. M(t) observada vs Hawkes — 3 CPs
# =============================================================================
plot_acumulada_hawkes_3cp <- function(n_modelo,
                                      samples_hawkes = samples_hawkes3,
                                      di,
                                      y_name     = "M(t)",
                                      col_obs    = "black",
                                      col_hawkes = "blue",
                                      col_banda  = "lightblue",
                                      col_cp     = "purple") {
  d  <- di[[n_modelo]]
  su <- samples_hawkes[[n_modelo]]$BUGSoutput$summary
  sl <- samples_hawkes[[n_modelo]]$BUGSoutput$sims.list
  alpha <- su[c("alpha[1]","alpha[2]","alpha[3]","alpha[4]"), "mean"]
  beta  <- su[c("beta[1]", "beta[2]", "beta[3]", "beta[4]"),  "mean"]
  n_par <- su["n",     "mean"]
  delta <- su["delta", "mean"]
  tau   <- su[c("tau[1]","tau[2]","tau[3]"), "mean"]
  t_grid   <- 1:1096
  acum_obs <- sapply(t_grid, function(i) sum(d <= i))
  cat(sprintf("Calculando bandas Hawkes 3CP modelo %d...\n", n_modelo))
  banda <- calc_M_hawkes_banda_3cp(
    t_grid, d, tau_sims=sl$tau,
    alpha_sims=sl$alpha, beta_sims=sl$beta,
    n_sims=sl$n, delta_sims=sl$delta)
  acum_est <- banda$mean
  plot(t_grid, acum_obs, type="l", lwd=2, lty="dashed", col=col_obs,
       xlab="Days", ylab=y_name,
       main=paste("Cumulative M(t) Hawkes 3CP:", y_name),
       ylim=c(0, max(c(acum_obs,banda$hi), na.rm=TRUE)*1.05))
  polygon(c(t_grid,rev(t_grid)), c(banda$hi,rev(banda$lo)),
          col=adjustcolor(col_banda, alpha.f=0.4), border=NA)
  lines(t_grid, banda$mean, lwd=2, col=col_hawkes)
  lines(t_grid, banda$lo, lwd=1, col=col_hawkes, lty="dashed")
  lines(t_grid, banda$hi, lwd=1, col=col_hawkes, lty="dashed")
  abline(v=tau, col=col_cp, lty="dotted", lwd=1.5)
  legend("topleft", legend=c("Observed","Hawkes M(t)","95% CI","CPs tau1/2/3"),
         col=c(col_obs,col_hawkes,col_hawkes,col_cp),
         lty=c("dashed","solid","dashed","dotted"), lwd=2, cex=0.8)
  sdm <- sum(abs(acum_obs[d] - banda$mean[d]))
  cat(sprintf("SDM modelo %d: %.4f\n", n_modelo, sdm))
  return(invisible(sdm))
}

# =============================================================================
# 4. M(t) observada vs NHPP vs Hawkes — 3 CPs
# =============================================================================
plot_acumulada_comparacion_3cp <- function(n_modelo,
                                           samples_nhpp   = samples_NHPP3,
                                           samples_hawkes = samples_hawkes3,
                                           di,
                                           y_name           = "M(t)",
                                           col_obs          = "black",
                                           col_nhpp         = "blue",
                                           col_hawkes       = "red",
                                           col_banda_nhpp   = "lightblue",
                                           col_banda_hawkes = "lightsalmon",
                                           col_cp           = "purple") {
  d      <- di[[n_modelo]]
  t_grid <- 1:1096
  acum_obs <- sapply(t_grid, function(i) sum(d <= i))
  # NHPP 3CP
  su_nhpp <- samples_nhpp[[n_modelo]]$BUGSoutput$summary
  sl_nhpp <- samples_nhpp[[n_modelo]]$BUGSoutput$sims.list
  a_n   <- su_nhpp[c("alpha[1]","alpha[2]","alpha[3]","alpha[4]"), "mean"]
  b_n   <- su_nhpp[c("beta[1]", "beta[2]", "beta[3]", "beta[4]"),  "mean"]
  tau_n <- su_nhpp[c("tau[1]","tau[2]","tau[3]"), "mean"]
  cat(sprintf("Calculando bandas NHPP 3CP modelo %d...\n", n_modelo))
  banda_nhpp <- calc_M_nhpp_banda_3cp(t_grid, sl_nhpp$tau, sl_nhpp$alpha, sl_nhpp$beta)
  acum_nhpp <- banda_nhpp$mean
  # Hawkes 3CP
  su_hwk <- samples_hawkes[[n_modelo]]$BUGSoutput$summary
  sl_hwk <- samples_hawkes[[n_modelo]]$BUGSoutput$sims.list
  a_h   <- su_hwk[c("alpha[1]","alpha[2]","alpha[3]","alpha[4]"), "mean"]
  b_h   <- su_hwk[c("beta[1]", "beta[2]", "beta[3]", "beta[4]"),  "mean"]
  n_h   <- su_hwk["n",     "mean"]
  d_h   <- su_hwk["delta", "mean"]
  tau_h <- su_hwk[c("tau[1]","tau[2]","tau[3]"), "mean"]
  cat(sprintf("Calculando bandas Hawkes 3CP modelo %d...\n", n_modelo))
  banda_hwk <- calc_M_hawkes_banda_3cp(
    t_grid, d, sl_hwk$tau,
    sl_hwk$alpha, sl_hwk$beta, sl_hwk$n, sl_hwk$delta)
  acum_hwk <- banda_hwk$mean
  ylim_max <- max(c(acum_obs,banda_nhpp$hi,banda_hwk$hi), na.rm=TRUE)*1.05
  plot(t_grid, acum_obs, type="l", lwd=2, lty="dashed", col=col_obs,
       xlab="Days", ylab=y_name,
       main=paste("M(t) comparacion 3CP:", y_name), ylim=c(0, ylim_max))
  polygon(c(t_grid,rev(t_grid)), c(banda_nhpp$hi,rev(banda_nhpp$lo)),
          col=adjustcolor(col_banda_nhpp, alpha.f=0.4), border=NA)
  lines(t_grid, banda_nhpp$mean, lwd=2, col=col_nhpp)
  lines(t_grid, banda_nhpp$lo, lwd=1, col=col_nhpp, lty="dashed")
  lines(t_grid, banda_nhpp$hi, lwd=1, col=col_nhpp, lty="dashed")
  polygon(c(t_grid,rev(t_grid)), c(banda_hwk$hi,rev(banda_hwk$lo)),
          col=adjustcolor(col_banda_hawkes, alpha.f=0.4), border=NA)
  lines(t_grid, banda_hwk$mean, lwd=2, col=col_hawkes)
  lines(t_grid, banda_hwk$lo, lwd=1, col=col_hawkes, lty="dashed")
  lines(t_grid, banda_hwk$hi, lwd=1, col=col_hawkes, lty="dashed")
  abline(v=tau_h, col=col_cp, lty="dotted", lwd=1.5)
  legend("topleft",
         legend=c("Observed","NHPP","NHPP 95% CI","Hawkes M4","Hawkes 95% CI","CPs tau1/2/3"),
         col=c(col_obs,col_nhpp,col_nhpp,col_hawkes,col_hawkes,col_cp),
         lty=c("dashed","solid","dashed","solid","dashed","dotted"), lwd=2, cex=0.8)
  sdm_nhpp   <- sum(abs(acum_obs[d] - banda_nhpp$mean[d]))
  sdm_hawkes <- sum(abs(acum_obs[d] - banda_hwk$mean[d]))
  cat(sprintf("Modelo %d — SDM NHPP: %.4f  |  SDM Hawkes: %.4f\n",
              n_modelo, sdm_nhpp, sdm_hawkes))
  return(invisible(list(sdm_nhpp=sdm_nhpp, sdm_hawkes=sdm_hawkes)))
}

# =============================================================================
# USO
# =============================================================================
load("NHPP_3cp.RData"); load("hawkes_3cp.RData")
dir_out <- "08_resultados_graficas/nhpp_vs_hawkes_3cp"
y_names <- c("PM2.5 Col 37 Bog","PM2.5 Col 37 Med","PM2.5 WHO 25 Bog","PM2.5 WHO 25 Med")

png(file.path(dir_out,"lambda_hawkes_3cp.png"), width=2400, height=1800, res=200)
par(mfrow=c(2,2), mar=c(4,4,2,1))
for (i in 1:4) plot_lambda_hawkes_3cp(n_modelo=i, di=di, y_name=y_names[i], grilla_paso=5)
par(mfrow=c(1,1)); dev.off()

png(file.path(dir_out,"baseline_nhpp_vs_hawkes_3cp.png"), width=2400, height=1800, res=200)
par(mfrow=c(2,2), mar=c(4,4,2,1))
for (i in 1:4) plot_baseline_comparison_3cp(n_modelo=i, y_name=y_names[i])
par(mfrow=c(1,1)); dev.off()

png(file.path(dir_out,"acumulada_hawkes_3cp.png"), width=2400, height=1800, res=200)
par(mfrow=c(2,2), mar=c(4,4,2,1))
for (i in 1:4) plot_acumulada_hawkes_3cp(n_modelo=i, di=di, y_name=y_names[i])
par(mfrow=c(1,1)); dev.off()

png(file.path(dir_out,"comparacion_Mt_nhpp_vs_hawkes_3cp.png"), width=2400, height=1800, res=200)
par(mfrow=c(2,2), mar=c(4,4,2,1))
for (i in 1:4) plot_acumulada_comparacion_3cp(n_modelo=i, di=di, y_name=y_names[i])
par(mfrow=c(1,1)); dev.off()

cat("\nGraficas guardadas en:", dir_out, "\n")
