# Visualización Hawkes M1 - 0 Change Points

calc_lambda_hawkes <- function(t_grid, d, alpha, beta, n, delta) {
  eta <- n * delta
  sapply(t_grid, function(ti) {
    mu_t <- (alpha / beta) * (ti / beta)^(alpha - 1)
    prev <- d[d < ti]
    kern <- if (length(prev) > 0) eta * sum(exp(-delta * (ti - prev))) else 0
    mu_t + kern
  })
}

# -----------------------------------------------------------------------------
# Auxiliar: M(t) Hawkes exacta — media posterior puntual
# M(t) = (t/beta)^alpha + (eta/delta) * sum_{d_k<t} (1 - exp(-delta*(t-d_k)))
# -----------------------------------------------------------------------------
calc_M_hawkes_exacta <- function(t_eval, d, alpha, beta, n, delta) {
  eta <- n * delta
  sapply(t_eval, function(tt) {
    prev <- d[d < tt]
    (tt / beta)^alpha +
      if (length(prev) > 0) (eta / delta) * sum(1 - exp(-delta * (tt - prev))) else 0
  })
}

# -----------------------------------------------------------------------------
# Auxiliar: bandas de credibilidad M(t) Hawkes — propaga incertidumbre MCMC
# Devuelve lista con mean, lo (2.5%), hi (97.5%)
# -----------------------------------------------------------------------------
calc_M_hawkes_banda <- function(t_eval, d, alpha_sims, beta_sims,
                                n_sims, delta_sims) {
  n_muestras <- length(n_sims)
  m_matrix   <- matrix(NA, nrow = n_muestras, ncol = length(t_eval))
  
  for (s in 1:n_muestras) {
    a_s     <- if (is.matrix(alpha_sims)) alpha_sims[s, 1] else alpha_sims[s]
    b_s     <- if (is.matrix(beta_sims))  beta_sims[s, 1]  else beta_sims[s]
    eta_s   <- n_sims[s] * delta_sims[s]
    delta_s <- delta_sims[s]
    
    m_matrix[s, ] <- sapply(t_eval, function(tt) {
      prev <- d[d < tt]
      (tt / b_s)^a_s +
        if (length(prev) > 0) {
          (eta_s / delta_s) * sum(1 - exp(-delta_s * (tt - prev)))
        } else { 0 }
    })
  }
  
  list(
    mean = apply(m_matrix, 2, mean),
    lo   = apply(m_matrix, 2, quantile, 0.025),
    hi   = apply(m_matrix, 2, quantile, 0.975)
  )
}

# -----------------------------------------------------------------------------
# Auxiliar: bandas de credibilidad M(t) NHPP — propaga incertidumbre MCMC
# M(t) = (t/beta)^alpha — forma cerrada, mas rapido que Hawkes
# Devuelve lista con mean, lo (2.5%), hi (97.5%)
# -----------------------------------------------------------------------------
calc_M_nhpp_banda <- function(t_eval, alpha_sims, beta_sims) {
  n_muestras <- length(alpha_sims)
  m_matrix   <- matrix(NA, nrow = n_muestras, ncol = length(t_eval))
  
  for (s in 1:n_muestras) {
    a_s <- if (is.matrix(alpha_sims)) alpha_sims[s, 1] else alpha_sims[s]
    b_s <- if (is.matrix(beta_sims))  beta_sims[s, 1]  else beta_sims[s]
    m_matrix[s, ] <- (t_eval / b_s)^a_s
  }
  
  list(
    mean = apply(m_matrix, 2, mean),
    lo   = apply(m_matrix, 2, quantile, 0.025),
    hi   = apply(m_matrix, 2, quantile, 0.975)
  )
}

# -----------------------------------------------------------------------------
# Auxiliar: baseline Weibull puntual
# -----------------------------------------------------------------------------
calc_mu_weibull <- function(t_grid, alpha, beta) {
  (alpha / beta) * (t_grid / beta)^(alpha - 1)
}

# =============================================================================
# 1. Intensidad condicional completa lambda(t | H_t)
# =============================================================================
plot_lambda_hawkes <- function(n_modelo,
                               samples_hawkes  = samples_hawkes0,
                               di,
                               y_name          = "Lambda",
                               grilla_paso     = 5,
                               mostrar_eventos = TRUE,
                               col_lambda      = "black",
                               col_mu          = "gray50",
                               col_eventos     = "red",
                               lwd_lambda      = 1.5) {
  
  d  <- di[[n_modelo]]
  su <- samples_hawkes[[n_modelo]]$BUGSoutput$summary
  
  alpha <- su["alpha", "mean"]
  beta  <- su["beta",  "mean"]
  n_par <- su["n",     "mean"]
  delta <- su["delta", "mean"]
  
  cat(sprintf("\nModelo %d: %s\n", n_modelo, y_name))
  cat(sprintf("  alpha=%.4f  beta=%.4f  n=%.4f  delta=%.4f  eta=%.4f\n",
              alpha, beta, n_par, delta, n_par * delta))
  cat(sprintf("  Vida media del kernel: %.2f dias\n", log(2) / delta))
  
  t_grid      <- seq(1, 1096, by = grilla_paso)
  lambda_full <- calc_lambda_hawkes(t_grid, d, alpha, beta, n_par, delta)
  mu_base     <- calc_mu_weibull(t_grid, alpha, beta)
  
  plot(t_grid, lambda_full,
       type = "l", lwd = lwd_lambda, col = col_lambda,
       xlab = "Days", ylab = y_name,
       main = paste("Hawkes M1-0CP:", y_name),
       ylim = c(0, max(lambda_full, na.rm = TRUE) * 1.1))
  lines(t_grid, mu_base, lwd = 1, col = col_mu, lty = "dashed")
  if (mostrar_eventos)
    rug(d, col = adjustcolor(col_eventos, alpha.f = 0.4), ticksize = 0.03)
  legend("topright",
         legend = c(expression(lambda^"*"*(t*"|"*H[t])),
                    expression(mu(t)~"baseline")),
         col    = c(col_lambda, col_mu),
         lty    = c("solid", "dashed"),
         lwd    = c(lwd_lambda, 1), cex = 0.8)
}

# =============================================================================
# 2. Comparacion baseline mu(t): NHPP vs Hawkes
# =============================================================================
plot_baseline_comparison <- function(n_modelo,
                                     samples_nhpp   = samples_resultados0,
                                     samples_hawkes = samples_hawkes0,
                                     y_name         = "mu(t)",
                                     col_nhpp       = "blue",
                                     col_hawkes     = "darkgreen") {
  
  t_grid  <- 1:1096
  su_nhpp <- samples_nhpp[[n_modelo]]$BUGSoutput$summary
  su_hwk  <- samples_hawkes[[n_modelo]]$BUGSoutput$summary
  
  mu_nhpp <- calc_mu_weibull(t_grid, su_nhpp["alpha","mean"], su_nhpp["beta","mean"])
  mu_hwk  <- calc_mu_weibull(t_grid, su_hwk["alpha","mean"],  su_hwk["beta","mean"])
  
  ylim_max <- max(c(mu_nhpp, mu_hwk), na.rm = TRUE) * 1.1
  
  plot(t_grid, mu_nhpp,
       type = "l", lwd = 2, col = col_nhpp,
       xlab = "Days", ylab = y_name,
       main = paste("Baseline NHPP vs Hawkes:", y_name),
       ylim = c(0, ylim_max))
  lines(t_grid, mu_hwk, lwd = 2, col = col_hawkes, lty = "dashed")
  legend("topright",
         legend = c("NHPP", "Hawkes M1"),
         col    = c(col_nhpp, col_hawkes),
         lty    = c("solid", "dashed"),
         lwd    = 2, cex = 0.8)
}

# =============================================================================
# 3. M(t) observada vs Hawkes con bandas de credibilidad
# =============================================================================
plot_acumulada_hawkes <- function(n_modelo,
                                  samples_hawkes = samples_hawkes0,
                                  di,
                                  y_name         = "M(t)",
                                  col_obs        = "black",
                                  col_hawkes     = "blue",
                                  col_banda      = "lightblue") {
  
  d  <- di[[n_modelo]]
  su <- samples_hawkes[[n_modelo]]$BUGSoutput$summary
  sl <- samples_hawkes[[n_modelo]]$BUGSoutput$sims.list
  
  alpha <- su["alpha", "mean"]
  beta  <- su["beta",  "mean"]
  n_par <- su["n",     "mean"]
  delta <- su["delta", "mean"]
  
  t_grid   <- 1:1096
  acum_obs <- sapply(t_grid, function(i) sum(d <= i))
  
  # Resumen posterior de la trayectoria M(t)
  cat(sprintf("Calculando bandas Hawkes modelo %d...\n", n_modelo))
  banda <- calc_M_hawkes_banda(t_grid, d,
                               sl$alpha, sl$beta, sl$n, sl$delta)
  acum_est <- banda$mean
  
  ylim_max <- max(c(acum_obs, banda$hi), na.rm = TRUE) * 1.05
  
  plot(t_grid, acum_obs,
       type = "l", lwd = 2, lty = "dashed", col = col_obs,
       xlab = "Days", ylab = y_name,
       main = paste("Cumulative M(t) Hawkes:", y_name),
       ylim = c(0, ylim_max))
  polygon(c(t_grid, rev(t_grid)),
          c(banda$hi, rev(banda$lo)),
          col = adjustcolor(col_banda, alpha.f = 0.4), border = NA)
  lines(t_grid, banda$mean, lwd = 2, col = col_hawkes)
  lines(t_grid, banda$lo, lwd = 1, col = col_hawkes, lty = "dashed")
  lines(t_grid, banda$hi, lwd = 1, col = col_hawkes, lty = "dashed")
  legend("topleft",
         legend = c("Observed", "Hawkes M(t)", "95% CI"),
         col    = c(col_obs, col_hawkes, col_hawkes),
         lty    = c("dashed", "solid", "dashed"),
         lwd    = 2, cex = 0.8)
  
  sdm <- sum(abs(acum_obs[d] - banda$mean[d]))
  cat(sprintf("SDM modelo %d: %.4f\n", n_modelo, sdm))
  return(invisible(sdm))
}

# =============================================================================
# 4. M(t) observada vs NHPP vs Hawkes — ambos con bandas de credibilidad
# =============================================================================
plot_acumulada_comparacion <- function(n_modelo,
                                       samples_nhpp   = samples_resultados0,
                                       samples_hawkes = samples_hawkes0,
                                       di,
                                       y_name         = "M(t)",
                                       col_obs        = "black",
                                       col_nhpp       = "blue",
                                       col_hawkes     = "red",
                                       col_banda_nhpp   = "lightblue",
                                       col_banda_hawkes = "lightsalmon") {
  
  d      <- di[[n_modelo]]
  t_grid <- 1:1096
  
  acum_obs <- sapply(t_grid, function(i) sum(d <= i))
  
  # --- NHPP ---
  su_nhpp <- samples_nhpp[[n_modelo]]$BUGSoutput$summary
  sl_nhpp <- samples_nhpp[[n_modelo]]$BUGSoutput$sims.list
  
  usar_summary_m <- TRUE
  banda_nhpp <- tryCatch(
    extract_M_nhpp_from_summary(su_nhpp),
    error = function(e) {
      usar_summary_m <<- FALSE
      NULL
    }
  )
  
  if (usar_summary_m) {
    cat(sprintf("Usando m del summary para NHPP modelo %d.\n", n_modelo))
  } else {
    cat(sprintf("Calculando bandas NHPP modelo %d desde sims.list...\n", n_modelo))
    banda_nhpp <- calc_M_nhpp_banda(t_grid, sl_nhpp$alpha, sl_nhpp$beta)
  }
  acum_nhpp <- banda_nhpp$mean
  
  # --- Hawkes ---
  su_hwk <- samples_hawkes[[n_modelo]]$BUGSoutput$summary
  sl_hwk <- samples_hawkes[[n_modelo]]$BUGSoutput$sims.list
  
  cat(sprintf("Calculando bandas Hawkes modelo %d...\n", n_modelo))
  banda_hwk <- calc_M_hawkes_banda(t_grid, d,
                                   sl_hwk$alpha, sl_hwk$beta,
                                   sl_hwk$n,     sl_hwk$delta)
  acum_hwk <- banda_hwk$mean
  
  ylim_max <- max(c(acum_obs, banda_nhpp$hi, banda_hwk$hi), na.rm = TRUE) * 1.05
  
  plot(t_grid, acum_obs,
       type = "l", lwd = 2, lty = "dashed", col = col_obs,
       xlab = "Days", ylab = y_name,
       main = paste("M(t) comparacion:", y_name),
       ylim = c(0, ylim_max))
  
  # Banda NHPP
  polygon(c(t_grid, rev(t_grid)),
          c(banda_nhpp$hi, rev(banda_nhpp$lo)),
          col = adjustcolor(col_banda_nhpp, alpha.f = 0.4), border = NA)
  lines(t_grid, banda_nhpp$mean, lwd = 2, col = col_nhpp)
  lines(t_grid, banda_nhpp$lo, lwd = 1, col = col_nhpp, lty = "dashed")
  lines(t_grid, banda_nhpp$hi, lwd = 1, col = col_nhpp, lty = "dashed")
  
  # Banda Hawkes
  polygon(c(t_grid, rev(t_grid)),
          c(banda_hwk$hi, rev(banda_hwk$lo)),
          col = adjustcolor(col_banda_hawkes, alpha.f = 0.4), border = NA)
  lines(t_grid, banda_hwk$mean, lwd = 2, col = col_hawkes)
  lines(t_grid, banda_hwk$lo, lwd = 1, col = col_hawkes, lty = "dashed")
  lines(t_grid, banda_hwk$hi, lwd = 1, col = col_hawkes, lty = "dashed")
  
  legend("topleft",
         legend = c("Observed", "NHPP", "NHPP 95% CI",
                    "Hawkes M1", "Hawkes 95% CI"),
         col    = c(col_obs, col_nhpp, col_nhpp, col_hawkes, col_hawkes),
         lty    = c("dashed", "solid", "dashed", "solid", "dashed"),
         lwd    = 2, cex = 0.8)
  
  # SDM ambos modelos
  sdm_nhpp   <- sum(abs(acum_obs[d] - acum_nhpp[d]))
  sdm_hawkes <- sum(abs(acum_obs[d] - acum_hwk[d]))
  cat(sprintf("Modelo %d — SDM NHPP: %.4f  |  SDM Hawkes: %.4f\n",
              n_modelo, sdm_nhpp, sdm_hawkes))
  
  return(invisible(list(sdm_nhpp = sdm_nhpp, sdm_hawkes = sdm_hawkes)))
}

# =============================================================================
# USO
# =============================================================================

load("NHPP_0cp.RData")
load("hawkes_0cp.RData")

dir_out <- "08_resultados_graficas/nhpp_vs_hawkes_0cp"

y_names <- c("PM2.5 Col 37 Bog", "PM2.5 Col 37 Med",
             "PM2.5 WHO 25 Bog", "PM2.5 WHO 25 Med")

# Índices por umbral
idx_col <- 1:2   # Col 37
idx_who <- 3:4   # WHO 25

# -----------------------------------------------------------------------------
# 1. Intensidad condicional lambda(t | H_t) — par Col y par WHO
# -----------------------------------------------------------------------------
png(file.path(dir_out, "lambda_hawkes_0cp_col.png"),
    width = 2400, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (i in idx_col) plot_lambda_hawkes(n_modelo = i, di = di, y_name = y_names[i], grilla_paso = 5)
par(mfrow = c(1, 1)); dev.off()

png(file.path(dir_out, "lambda_hawkes_0cp_who.png"),
    width = 2400, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (i in idx_who) plot_lambda_hawkes(n_modelo = i, di = di, y_name = y_names[i], grilla_paso = 5)
par(mfrow = c(1, 1)); dev.off()

# -----------------------------------------------------------------------------
# 2. Comparacion baseline mu(t): NHPP vs Hawkes — par Col y par WHO
# -----------------------------------------------------------------------------
png(file.path(dir_out, "baseline_nhpp_vs_hawkes_0cp_col.png"),
    width = 2400, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (i in idx_col) plot_baseline_comparison(n_modelo = i, y_name = y_names[i])
par(mfrow = c(1, 1)); dev.off()

png(file.path(dir_out, "baseline_nhpp_vs_hawkes_0cp_who.png"),
    width = 2400, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (i in idx_who) plot_baseline_comparison(n_modelo = i, y_name = y_names[i])
par(mfrow = c(1, 1)); dev.off()

# -----------------------------------------------------------------------------
# 3. M(t) Hawkes con bandas de credibilidad — par Col y par WHO
# -----------------------------------------------------------------------------
png(file.path(dir_out, "acumulada_hawkes_0cp_col.png"),
    width = 2400, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (i in idx_col) plot_acumulada_hawkes(n_modelo = i, di = di, y_name = y_names[i])
par(mfrow = c(1, 1)); dev.off()

png(file.path(dir_out, "acumulada_hawkes_0cp_who.png"),
    width = 2400, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (i in idx_who) plot_acumulada_hawkes(n_modelo = i, di = di, y_name = y_names[i])
par(mfrow = c(1, 1)); dev.off()

# -----------------------------------------------------------------------------
# 4. Comparacion M(t) NHPP vs Hawkes — ambos con bandas — par Col y par WHO
# -----------------------------------------------------------------------------
png(file.path(dir_out, "comparacion_Mt_nhpp_vs_hawkes_0cp_col.png"),
    width = 2400, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (i in idx_col) plot_acumulada_comparacion(n_modelo = i, di = di, y_name = y_names[i])
par(mfrow = c(1, 1)); dev.off()

png(file.path(dir_out, "comparacion_Mt_nhpp_vs_hawkes_0cp_who.png"),
    width = 2400, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (i in idx_who) plot_acumulada_comparacion(n_modelo = i, di = di, y_name = y_names[i])
par(mfrow = c(1, 1)); dev.off()

cat("\nGraficas guardadas en:", dir_out, "\n")