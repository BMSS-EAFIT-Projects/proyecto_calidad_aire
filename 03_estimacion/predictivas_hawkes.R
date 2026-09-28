# Only edit the "CONFIGURATION" section.

library(ggplot2)
library(dplyr)
library(tidyr)


# =============================================================================
# CONFIGURATION
# =============================================================================

# Names of the .RData files
rdata_files <- c(
  "hawkes_0cp",
  "hawkes_1cp",
  "hawkes_2cp",
  "hawkes_3cp"
)

# Name of the R2jags object inside each .RData
obj_names <- rdata_files

# Vectors of observed exceedances
di <- list(excePM25col_Bog, excePM25col_Med, excePM25who_Bog, excePM25who_Med)
d_names <- c("PM2.5 Col37 Bog", "PM2.5 Col37 Med",
             "PM2.5 WHO Bog",   "PM2.5 WHO Med")

T_obs   <- 1096   # last observed day
H       <- 180    # forecast horizon
T_pred  <- T_obs + H

# Number of MCMC samples to use for the band
n_sims_use <- 4000

# Output directory for the charts
output_dir <- "08_resultados_graficas/predictiva_hawkes"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# HELPER FUNCTIONS — mean function M(t)
# =============================================================================

calc_M_hawkes <- function(t_eval, d_hist, alpha, beta, eta, delta) {
  sapply(t_eval, function(tt) {
    prev         <- d_hist[d_hist < tt]
    weibull_part <- (tt / beta)^alpha
    hawkes_part  <- if (length(prev) > 0)
      (eta / delta) * sum(1 - exp(-delta * (tt - prev)))
    else 0
    weibull_part + hawkes_part
  })
}

# M(t) for models with change-points
calc_M_hawkes_cp <- function(t_eval, d_hist, alpha_vec, beta_vec,
                             tau_vec, eta, delta) {
  n_seg <- length(alpha_vec)
  
  M_single <- function(tt) {
    if (n_seg == 1) {
      seg <- 1
    } else {
      seg <- min(findInterval(tt, c(0, tau_vec, Inf)), n_seg)
    }
    
    # Weibull part with continuity correction
    if (seg == 1) {
      w <- (tt / beta_vec[1])^alpha_vec[1]
    } else {
      m_tau1 <- (tau_vec[1] / beta_vec[1])^alpha_vec[1]
      if (seg == 2) {
        w <- m_tau1 +
          (tt         / beta_vec[2])^alpha_vec[2] -
          (tau_vec[1] / beta_vec[2])^alpha_vec[2]
      } else if (seg == 3) {
        delta2 <- (tau_vec[2] / beta_vec[2])^alpha_vec[2] -
          (tau_vec[1] / beta_vec[2])^alpha_vec[2]
        w <- m_tau1 + delta2 +
          (tt         / beta_vec[3])^alpha_vec[3] -
          (tau_vec[2] / beta_vec[3])^alpha_vec[3]
      } else {  # seg == 4 (3CP)
        delta2 <- (tau_vec[2] / beta_vec[2])^alpha_vec[2] -
          (tau_vec[1] / beta_vec[2])^alpha_vec[2]
        delta3 <- (tau_vec[3] / beta_vec[3])^alpha_vec[3] -
          (tau_vec[2] / beta_vec[3])^alpha_vec[3]
        w <- m_tau1 + delta2 + delta3 +
          (tt         / beta_vec[4])^alpha_vec[4] -
          (tau_vec[3] / beta_vec[4])^alpha_vec[4]
      }
    }
    
    # Hawkes part (memory of all observed events before tt)
    prev        <- d_hist[d_hist < tt]
    hawkes_part <- if (length(prev) > 0)
      (eta / delta) * sum(1 - exp(-delta * (tt - prev)))
    else 0
    
    w + hawkes_part
  }
  
  sapply(t_eval, M_single)
}


# =============================================================================
# PARAMETER EXTRACTION FROM sims.list
# =============================================================================

extraer_params <- function(sl, n_cp) {
  n_sims <- length(sl$eta)
  n_seg  <- n_cp + 1
  
  if (n_cp == 0) {
    alpha_mat <- matrix(as.numeric(sl$alpha), ncol = 1)
    beta_mat  <- matrix(as.numeric(sl$beta),  ncol = 1)
    tau_mat   <- NULL
  } else {
    alpha_mat <- if (is.matrix(sl$alpha)) sl$alpha else matrix(sl$alpha, ncol = n_seg)
    beta_mat  <- if (is.matrix(sl$beta))  sl$beta  else matrix(sl$beta,  ncol = n_seg)
    # tau: vector (1CP) or matrix (2CP, 3CP) — always coerce to matrix
    tau_mat   <- if (is.matrix(sl$tau)) sl$tau else matrix(sl$tau, ncol = n_cp)
  }
  
  list(
    n_sims    = n_sims,
    n_seg     = n_seg,
    alpha_mat = alpha_mat,
    beta_mat  = beta_mat,
    tau_mat   = tau_mat,
    eta_vec   = as.numeric(sl$eta),
    delta_vec = as.numeric(sl$delta)
  )
}


# =============================================================================
# OGATA THINNING
# =============================================================================

sim_hawkes_ogata <- function(T_obs, H, d_hist,
                             alpha_s, beta_s, eta_s, delta_s,
                             n_cp, tau_s) {
  
  # ---- Background rate λ₀(t) = d/dt [M_weibull(t)] for the last segment ----
  a_fc <- alpha_s[n_cp + 1]
  b_fc <- beta_s[n_cp + 1]
  weibull_rate <- function(t) (a_fc / b_fc) * (t / b_fc)^(a_fc - 1)
  
  # ---- Conditional intensity λ*(t | history) --------------------------------
  lambda_star <- function(t, history) {
    bg     <- weibull_rate(t)
    kernel <- if (length(history) > 0)
      eta_s * sum(exp(-delta_s * (t - history[history < t])))
    else 0
    bg + kernel
  }
  
  history    <- d_hist
  t_current  <- T_obs
  new_events <- numeric(0)
  T_end      <- T_obs + H
  
  while (TRUE) {
    wb_upper <- max(weibull_rate(t_current), weibull_rate(T_end))
    
    # Kernel part of the upper bound: non-increasing, so current value suffices
    kernel_upper <- if (length(history) > 0)
      eta_s * sum(exp(-delta_s * (t_current - history[history < t_current])))
    else 0
    
    lam_bar <- wb_upper + kernel_upper
    if (lam_bar <= 0) lam_bar <- 1e-9
    
    # Candidate arrival time
    t_cand <- t_current + rexp(1, rate = lam_bar)
    if (t_cand > T_end) break
    
    # True intensity at candidate — guaranteed ≤ lam_bar by construction
    lam_true <- lambda_star(t_cand, history)
    
    # Accept/reject (ratio always in [0, 1])
    if (runif(1) < lam_true / lam_bar) {
      history    <- c(history, t_cand)
      new_events <- c(new_events, t_cand)
    }
    
    t_current <- t_cand
  }
  
  new_events
}


# =============================================================================
# POSTERIOR PREDICTIVE BAND CALCULATION
# =============================================================================

calcular_bandas_pp <- function(d, params, t_grid_obs, t_grid_pred,
                               T_obs, H, n_sims_use) {
  n_cp  <- params$n_seg - 1
  n_all <- params$n_sims
  
  set.seed(42)
  idx <- sample(n_all, min(n_sims_use, n_all))
  
  n_obs  <- length(t_grid_obs)
  n_pred <- length(t_grid_pred)
  
  N_obs_mat  <- matrix(NA_real_, length(idx), n_obs)
  N_pred_mat <- matrix(NA_real_, length(idx), n_pred)
  M_obs_mat  <- matrix(NA_real_, length(idx), n_obs)
  M_pred_mat <- matrix(NA_real_, length(idx), n_pred)
  
  N_T_obs <- sum(d <= T_obs)
  
  for (k in seq_along(idx)) {
    s       <- idx[k]
    alpha_s <- as.numeric(params$alpha_mat[s, ])
    beta_s  <- as.numeric(params$beta_mat[s, ])
    eta_s   <- params$eta_vec[s]
    delta_s <- params$delta_vec[s]
    tau_s   <- if (!is.null(params$tau_mat)) as.numeric(params$tau_mat[s, ]) else numeric(0)
    
    # ---- In-sample: evaluate M(t | θ_s) and draw Poisson increments --------
    if (n_cp == 0) {
      M_obs <- calc_M_hawkes(t_grid_obs, d, alpha_s, beta_s, eta_s, delta_s)
    } else {
      M_obs <- calc_M_hawkes_cp(t_grid_obs, d, alpha_s, beta_s, tau_s, eta_s, delta_s)
    }
    M_obs_mat[k, ] <- M_obs
    
    delta_m_obs    <- pmax(diff(c(0, M_obs)), 0)
    N_obs_mat[k, ] <- cumsum(rpois(n_obs, lambda = delta_m_obs))
    
    # ---- Forecast: Ogata thinning 
    new_ev <- sim_hawkes_ogata(
      T_obs   = T_obs,
      H       = H,
      d_hist  = d,        
      alpha_s = alpha_s,
      beta_s  = beta_s,
      eta_s   = eta_s,
      delta_s = delta_s,
      n_cp    = n_cp,
      tau_s   = tau_s
    )
    
    # Count generated events at each forecast grid point
    N_pred_mat[k, ] <- N_T_obs + sapply(t_grid_pred, function(t) sum(new_ev <= t))
    
    d_ext <- c(d, new_ev)
    if (n_cp == 0) {
      M_pred <- calc_M_hawkes(t_grid_pred, d_ext, alpha_s, beta_s, eta_s, delta_s)
    } else {
      M_pred <- calc_M_hawkes_cp(t_grid_pred, d_ext, alpha_s, beta_s, tau_s, eta_s, delta_s)
    }
    M_pred_mat[k, ] <- M_pred
  }
  
  list(
    obs_pp  = list(
      mean  = colMeans(N_obs_mat,  na.rm = TRUE),
      lower = apply(N_obs_mat,  2, quantile, 0.025, na.rm = TRUE),
      upper = apply(N_obs_mat,  2, quantile, 0.975, na.rm = TRUE)
    ),
    pred_pp = list(
      mean  = colMeans(N_pred_mat,  na.rm = TRUE),
      lower = apply(N_pred_mat,  2, quantile, 0.025, na.rm = TRUE),
      upper = apply(N_pred_mat,  2, quantile, 0.975, na.rm = TRUE)
    ),
    obs_M  = list(
      mean  = colMeans(M_obs_mat,  na.rm = TRUE),
      lower = apply(M_obs_mat,  2, quantile, 0.025, na.rm = TRUE),
      upper = apply(M_obs_mat,  2, quantile, 0.975, na.rm = TRUE)
    ),
    pred_M = list(
      mean  = colMeans(M_pred_mat,  na.rm = TRUE),
      lower = apply(M_pred_mat,  2, quantile, 0.025, na.rm = TRUE),
      upper = apply(M_pred_mat,  2, quantile, 0.975, na.rm = TRUE)
    )
  )
}


# =============================================================================
# CHART FUNCTION  (posterior predictive version)
# =============================================================================

graficar_predictiva_pp <- function(d, bandas,
                                   t_grid_obs, t_grid_pred,
                                   titulo,
                                   T_obs, output_file) {
  
  # Observed cumulative count
  acum_obs    <- sapply(1:T_obs, function(i) sum(d <= i))
  df_obs_real <- data.frame(t = 1:T_obs, M = acum_obs)
  
  # --- Posterior predictive bands for N(t) ---
  df_pp_obs <- data.frame(
    t     = t_grid_obs,
    mean  = bandas$obs_pp$mean,
    lower = bandas$obs_pp$lower,
    upper = bandas$obs_pp$upper
  )
  df_pp_pred <- data.frame(
    t     = t_grid_pred,
    mean  = bandas$pred_pp$mean,
    lower = bandas$pred_pp$lower,
    upper = bandas$pred_pp$upper
  )
  
  # --- Posterior mean of M(t | θ) for reference ---
  df_M_obs  <- data.frame(t = t_grid_obs,  mean = bandas$obs_M$mean)
  df_M_pred <- data.frame(t = t_grid_pred, mean = bandas$pred_M$mean)
  
  p <- ggplot() +
    
    # ---- In-sample band and mean: no color aes → excluded from legend --------
  geom_ribbon(data = df_pp_obs,
              aes(x = t, ymin = lower, ymax = upper),
              fill = "#3182bd", alpha = 0.18, show.legend = FALSE) +
    geom_line(data = df_pp_obs,
              aes(x = t, y = mean),
              color = "#3182bd", linewidth = 0.85, show.legend = FALSE) +
    # ---- Forecast band and mean: in legend -----------------------------------
  geom_ribbon(data = df_pp_pred,
              aes(x = t, ymin = lower, ymax = upper),
              fill = "#e6550d", alpha = 0.18) +
    geom_line(data = df_pp_pred,
              aes(x = t, y = mean, color = "Post. pred. mean N(t) — forecast"),
              linewidth = 0.85, linetype = "dashed") +
    
    # ---- Observed cumulative count: in legend --------------------------------
  geom_line(data = df_obs_real,
            aes(x = t, y = M, color = "Observed N(t)"),
            linewidth = 0.7) +
    
    # ---- Separator ----------------------------------------------------------
  geom_vline(xintercept = T_obs,
             linetype = "dotted", color = "grey40", linewidth = 0.6) +
    annotate("text", x = T_obs + 5,
             y = min(df_pp_obs$lower, na.rm = TRUE),
             label = "t = T", hjust = 0, size = 3, color = "grey40") +
    
    scale_color_manual(
      name   = NULL,
      values = c(
        "Observed N(t)"                      = "black",
        "Post. pred. mean N(t) — forecast"   = "#e6550d"
      )
    ) +
    labs(
      title    = titulo,
      x        = "Days",
      y        = "Cumulative exceedance count"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      legend.key.width = unit(1.5, "cm"),
      plot.title       = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(nrow = 1))
  
  ggsave(output_file, plot = p, width = 8, height = 5.5, dpi = 300)
  message("  Saved: ", output_file)
  invisible(p)
}


# =============================================================================
# MAIN LOOP: 4 datasets × 4 CP models
# =============================================================================

# Time grids
t_grid_obs  <- seq(1, T_obs,  by = 2)
t_grid_pred <- seq(T_obs, T_pred, by = 1)

cp_labels <- c("0CP", "1CP", "2CP", "3CP")

for (cp_idx in 1:4) {
  n_cp       <- cp_idx - 1
  rdata_path <- paste0(rdata_files[cp_idx], ".RData")
  
  message("\n=== Loading: ", rdata_files[cp_idx], ".RData ===")
  env_tmp <- new.env()
  load(rdata_path, envir = env_tmp)
  
  # Try with the expected name; if not found, use the first object in env
  obj_name <- obj_names[cp_idx]
  if (!exists(obj_name, envir = env_tmp)) {
    obj_name <- ls(env_tmp)[1]
    message("  Object auto-detected: ", obj_name)
  }
  hawkes_obj <- get(obj_name, envir = env_tmp)
  
  # List of 4 models
  if (!is.null(hawkes_obj$BUGSoutput)) {
    hawkes_lista <- list(hawkes_obj)
    n_mods <- 1
    message("  Structure: single object (1 dataset)")
  } else if (is.list(hawkes_obj) && length(hawkes_obj) >= 1 &&
             !is.null(hawkes_obj[[1]]$BUGSoutput)) {
    hawkes_lista <- hawkes_obj
    n_mods <- length(hawkes_obj)
    message("  Structure: list of ", n_mods, " datasets")
  } else {
    warning("Unrecognised structure for ", rdata_files[cp_idx], " — skipping.")
    next
  }
  
  for (mod_idx in seq_len(n_mods)) {
    d        <- di[[mod_idx]]
    mod_name <- d_names[mod_idx]
    slug     <- gsub("[^A-Za-z0-9]", "_", mod_name)
    
    message("  Dataset: ", mod_name)
    
    sl     <- hawkes_lista[[mod_idx]]$BUGSoutput$sims.list
    params <- extraer_params(sl, n_cp)
    
    # Compute full posterior predictive bands (Ogata thinning for forecast)
    br   <- round(mean(params$eta_vec / params$delta_vec), 3)
    hl   <- round(mean(log(2) / params$delta_vec), 1)
    message("    Computing posterior predictive bands (", n_sims_use,
            " samples) — Ogata thinning ...")
    message("    Branching ratio eta/delta = ", br,
            "  |  Kernel half-life = ", hl, " days")
    bandas <- calcular_bandas_pp(d, params, t_grid_obs, t_grid_pred,
                                 T_obs = T_obs, H = H,
                                 n_sims_use = n_sims_use)
    
    # Generate chart
    output_file <- file.path(
      output_dir,
      sprintf("postpred_hawkes_%s_%s.png", cp_labels[cp_idx], slug)
    )
    
    graficar_predictiva_pp(
      d           = d,
      bandas      = bandas,
      t_grid_obs  = t_grid_obs,
      t_grid_pred = t_grid_pred,
      titulo      = sprintf("Posterior Predictive — Hawkes %s | %s", cp_labels[cp_idx], mod_name),
      T_obs       = T_obs,
      output_file = output_file
    )
  }
  
  rm(env_tmp)
  gc()
}

message("\n=== Process complete. Charts saved in: ", output_dir, " ===")