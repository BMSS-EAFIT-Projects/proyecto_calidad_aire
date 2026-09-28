library(loo)

load("hawkes_0cp.RData")
load("hawkes_1cp.RData")
load("hawkes_2cp.RData")
load("hawkes_3cp.RData")

`%+%` <- paste0

d_labels_tex <- c(
  "PM$_{2.5}$ Colombia 37 --- Bogot\\'{a}",
  "PM$_{2.5}$ Colombia 37 --- Medell\\'{i}n",
  "PM$_{2.5}$ WHO 25 --- Bogot\\'{a}",
  "PM$_{2.5}$ WHO 25 --- Medell\\'{i}n"
)

di <- list(excePM25col_Bog, excePM25col_Med,
           excePM25who_Bog, excePM25who_Med)

t_total      <- 1096
samples_list <- list(samples_hawkes0, samples_hawkes1,
                     samples_hawkes2, samples_hawkes3)
n_seg        <- c(1, 2, 3, 4)

# Priors Hawkes
prior_alpha <- "$\\text{Unif}(0,\\,3)$"
prior_beta  <- "$\\text{Unif}(0,\\,400)$"
prior_tau   <- c("$\\text{Unif}(0,\\,365)$",
                 "$\\text{Unif}(365,\\,730)$",
                 "$\\text{Unif}(730,\\,1096)$")
prior_n     <- "$\\text{Unif}(0,\\,1)$"
prior_delta <- "$\\text{Unif}(0,\\,1)$"

# Carpetas
base_dir <- "10_resultados_tabulares"
cp_dirs  <- base_dir %+% "/resultados_" %+% 0:3 %+% "cp"
for (d in c(base_dir, cp_dirs)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

ncols_of <- function(col_def) nchar(gsub("[^lrc]", "", col_def))

write_tabular <- function(caption, label, col_def, header, body, path,
                          footnote = NULL) {
  nc <- ncols_of(col_def)
  fn_lines <- if (!is.null(footnote))
    c("\\\\[-0.5em]",
      "\\multicolumn{" %+% nc %+% "}{l}{\\footnotesize " %+% footnote %+% "}")
  else character(0)
  
  lines <- c(
    "% !TeX root = main.tex",
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{" %+% caption %+% "}",
    "\\label{" %+% label %+% "}",
    "\\begin{tabular}{" %+% col_def %+% "}",
    "\\toprule",
    header,
    "\\midrule",
    body,
    "\\bottomrule",
    fn_lines,
    "\\end{tabular}",
    "\\end{table}"
  )
  writeLines(lines, con = path)
  cat("Guardado:", path, "\n")
}

write_longtable <- function(caption, label, col_def, header, body, path) {
  nc <- ncols_of(col_def)
  lines <- c(
    "% !TeX root = main.tex",
    "% Requiere \\usepackage{booktabs,longtable} en el preambulo",
    "\\begin{longtable}{" %+% col_def %+% "}",
    "\\caption{" %+% caption %+% "} \\label{" %+% label %+% "} \\\\",
    "\\toprule",
    header,
    "\\midrule",
    "\\endfirsthead",
    "",
    "\\multicolumn{" %+% nc %+% "}{l}{\\tablename\\ \\thetable{} -- continuaci\\'{o}n} \\\\",
    "\\toprule",
    header,
    "\\midrule",
    "\\endhead",
    "",
    "\\midrule",
    "\\multicolumn{" %+% nc %+% "}{r}{Contin\\'{u}a en la siguiente p\\'{a}gina} \\\\",
    "\\endfoot",
    "",
    "\\bottomrule",
    "\\endlastfoot",
    "",
    body,
    "\\end{longtable}"
  )
  writeLines(lines, con = path)
  cat("Guardado:", path, "\n")
}

# =============================================================================
# HELPER: M(t) con puntos de cambio y correccion de continuidad
# Identica a calc_M_hawkes_cp del script de predictivas.
# Para 0CP: alpha_vec y beta_vec son vectores de longitud 1, tau_vec = numeric(0)
# =============================================================================

calc_M_hawkes_cp <- function(t_eval, d_hist, alpha_vec, beta_vec,
                             tau_vec, eta, delta) {
  n_seg_loc <- length(alpha_vec)
  
  M_single <- function(tt) {
    # --- Segmento al que pertenece tt ---
    if (n_seg_loc == 1) {
      seg <- 1
    } else {
      seg <- min(findInterval(tt, c(0, tau_vec, Inf)), n_seg_loc)
    }
    
    # --- Parte Weibull con correccion de continuidad ---
    if (seg == 1) {
      w <- (tt / beta_vec[1])^alpha_vec[1]
    } else if (seg == 2) {
      m_tau1 <- (tau_vec[1] / beta_vec[1])^alpha_vec[1]
      w <- m_tau1 +
        (tt         / beta_vec[2])^alpha_vec[2] -
        (tau_vec[1] / beta_vec[2])^alpha_vec[2]
    } else if (seg == 3) {
      m_tau1 <- (tau_vec[1] / beta_vec[1])^alpha_vec[1]
      d2     <- (tau_vec[2] / beta_vec[2])^alpha_vec[2] -
        (tau_vec[1] / beta_vec[2])^alpha_vec[2]
      w <- m_tau1 + d2 +
        (tt         / beta_vec[3])^alpha_vec[3] -
        (tau_vec[2] / beta_vec[3])^alpha_vec[3]
    } else {  # seg == 4  (3CP)
      m_tau1 <- (tau_vec[1] / beta_vec[1])^alpha_vec[1]
      d2     <- (tau_vec[2] / beta_vec[2])^alpha_vec[2] -
        (tau_vec[1] / beta_vec[2])^alpha_vec[2]
      d3     <- (tau_vec[3] / beta_vec[3])^alpha_vec[3] -
        (tau_vec[2] / beta_vec[3])^alpha_vec[3]
      w <- m_tau1 + d2 + d3 +
        (tt         / beta_vec[4])^alpha_vec[4] -
        (tau_vec[3] / beta_vec[4])^alpha_vec[4]
    }
    
    # --- Parte Hawkes (kernel exponencial) ---
    prev        <- d_hist[d_hist < tt]
    hawkes_part <- if (length(prev) > 0)
      (eta / delta) * sum(1 - exp(-delta * (tt - prev)))
    else 0
    
    w + hawkes_part
  }
  
  sapply(t_eval, M_single)
}

# =============================================================================
# 1. Tabla de parametros — longtable, convencion articulo base
# =============================================================================

for (cp in 1:4) {
  ncp  <- cp - 1
  body <- character(0)
  
  for (mod in 1:4) {
    su     <- samples_list[[cp]][[mod]]$BUGSoutput$summary
    nombre <- d_labels_tex[mod]
    first  <- TRUE
    
    emit <- function(model_cell, param_tex, prior_tex, row_su) {
      mn  <- sprintf("%.2f", row_su[1, "mean"])
      sdd <- sprintf("%.2f", row_su[1, "sd"])
      lo  <- sprintf("%.2f", row_su[1, "2.5%"])
      hi  <- sprintf("%.2f", row_su[1, "97.5%"])
      ci  <- "(" %+% lo %+% ", " %+% hi %+% ")"
      paste(model_cell, "&", param_tex, "&", prior_tex, "&",
            mn, "&", sdd, "&", ci, "\\\\")
    }
    
    for (i in 1:n_seg[cp]) {
      a_rows <- su[grep("^alpha\\[" %+% i %+% "\\]", rownames(su)), , drop = FALSE]
      b_rows <- su[grep("^beta\\["  %+% i %+% "\\]", rownames(su)), , drop = FALSE]
      if (nrow(a_rows) == 0) a_rows <- su[grep("^alpha$", rownames(su)), , drop = FALSE]
      if (nrow(b_rows) == 0) b_rows <- su[grep("^beta$",  rownames(su)), , drop = FALSE]
      
      mc <- if (first) nombre else ""
      body  <- c(body, emit(mc, "$\\alpha_{" %+% i %+% "}$", prior_alpha, a_rows))
      first <- FALSE
      body  <- c(body, emit("",  "$\\beta_{" %+%  i %+% "}$", prior_beta,  b_rows))
    }
    
    if (ncp > 0) {
      for (j in 1:ncp) {
        t_rows <- su[grep("^tau\\[" %+% j %+% "\\]", rownames(su)), , drop = FALSE]
        if (nrow(t_rows) == 0) t_rows <- su[grep("^tau$", rownames(su)), , drop = FALSE]
        body <- c(body, emit("", "$\\tau_{" %+% j %+% "}$", prior_tau[j], t_rows))
      }
    }
    
    for (par_info in list(
      list(pat = "^n$",     tex = "$n$",       prior = prior_n),
      list(pat = "^delta$", tex = "$\\delta$",  prior = prior_delta),
      list(pat = "^eta$",   tex = "$\\eta$",    prior = "derivado: $n \\cdot \\delta$"))) {
      p_rows <- su[grep(par_info$pat, rownames(su)), , drop = FALSE]
      if (nrow(p_rows) > 0)
        body <- c(body, emit("", par_info$tex, par_info$prior, p_rows))
    }
    
    body <- c(body, "\\midrule")
  }
  body <- body[-length(body)]
  
  col_def <- "llllrrr"
  header  <- paste("Model & Parameters & Prior Distribution &",
                   "Mean & SD & Credible intervals 95\\% \\\\")
  caption <- "Prior distributions and estimated parameters --- Hawkes M1 " %+% ncp %+% "CP"
  label   <- "tab:hawkes_params_" %+% ncp %+% "cp"
  path    <- cp_dirs[cp] %+% "/parametros_hawkes_" %+% ncp %+% "cp.tex"
  
  write_longtable(caption, label, col_def, header, body, path)
}

# =============================================================================
# 2. Tabla WAIC
# =============================================================================

WAIC_mat <- matrix(NA, 4, 4)
for (cp in 1:4)
  for (mod in 1:4) {
    ll <- samples_list[[cp]][[mod]]$BUGSoutput$sims.list$phi2
    WAIC_mat[mod, cp] <- round(
      waic(ll)$estimates["elpd_waic", "Estimate"] * -2, 2)
  }

{
  header <- "Dataset & 0\\,CP & 1\\,CP & 2\\,CP & 3\\,CP \\\\"
  body   <- sapply(1:4, function(mod)
    d_labels_tex[mod] %+% " & " %+%
      paste(sprintf("%.2f", WAIC_mat[mod, ]), collapse = " & ") %+% " \\\\")
  
  write_tabular(
    caption = "WAIC por modelo y n\\'{u}mero de puntos de cambio --- Hawkes M1",
    label   = "tab:hawkes_waic",
    col_def = "lrrrr",
    header  = header,
    body    = body,
    path    = base_dir %+% "/WAIC_hawkes_todos.tex"
  )
  
  for (cp in 1:4) {
    ncp <- cp - 1
    body_cp <- sapply(1:4, function(mod)
      d_labels_tex[mod] %+% " & " %+% sprintf("%.2f", WAIC_mat[mod, cp]) %+% " \\\\")
    write_tabular(
      caption = "WAIC --- Hawkes M1 " %+% ncp %+% "CP",
      label   = "tab:hawkes_waic_" %+% ncp %+% "cp",
      col_def = "lr",
      header  = "Dataset & WAIC (" %+% ncp %+% "\\,CP) \\\\",
      body    = body_cp,
      path    = cp_dirs[cp] %+% "/WAIC_hawkes_" %+% ncp %+% "cp.tex"
    )
  }
}

# =============================================================================
# 3. Tabla DIC
# =============================================================================

DIC_mat <- matrix(NA, 4, 4)
for (cp in 1:4)
  for (mod in 1:4)
    DIC_mat[mod, cp] <- round(
      samples_list[[cp]][[mod]]$BUGSoutput$DIC, 2)

{
  header <- "Dataset & 0\\,CP & 1\\,CP & 2\\,CP & 3\\,CP \\\\"
  body   <- sapply(1:4, function(mod)
    d_labels_tex[mod] %+% " & " %+%
      paste(sprintf("%.2f", DIC_mat[mod, ]), collapse = " & ") %+% " \\\\")
  
  write_tabular(
    caption = "DIC por modelo y n\\'{u}mero de puntos de cambio --- Hawkes M1",
    label   = "tab:hawkes_dic",
    col_def = "lrrrr",
    header  = header,
    body    = body,
    path    = base_dir %+% "/DIC_hawkes_todos.tex"
  )
  
  for (cp in 1:4) {
    ncp <- cp - 1
    body_cp <- sapply(1:4, function(mod)
      d_labels_tex[mod] %+% " & " %+% sprintf("%.2f", DIC_mat[mod, cp]) %+% " \\\\")
    write_tabular(
      caption = "DIC --- Hawkes M1 " %+% ncp %+% "CP",
      label   = "tab:hawkes_dic_" %+% ncp %+% "cp",
      col_def = "lr",
      header  = "Dataset & DIC (" %+% ncp %+% "\\,CP) \\\\",
      body    = body_cp,
      path    = cp_dirs[cp] %+% "/DIC_hawkes_" %+% ncp %+% "cp.tex"
    )
  }
    }

# =============================================================================
# 5. Tabla SDM — CORREGIDO: usa calc_M_hawkes_cp con correccion de continuidad
#    SDM = (1/K) * sum_i |N(t_i) - M_hat(t_i)|
#    donde M_hat(t_i) es la media posterior del compensador evaluada
#    en cada tiempo de evento, propagando correctamente los puntos de cambio.
# =============================================================================

calc_sdm_hawkes <- function(d, sl, n_cp) {
  n_sims    <- length(sl$n)
  n_seg_loc <- n_cp + 1
  K         <- length(d)
  
  # N(t_i) = rango del evento i-esimo = i (eventos ordenados)
  acum_eventos <- seq_len(K)
  
  # Matriz de M(t_i) para cada muestra MCMC
  m_matrix <- matrix(NA, nrow = n_sims, ncol = K)
  
  for (s in 1:n_sims) {
    # Extraer parametros de la muestra s
    if (n_seg_loc == 1) {
      alpha_s <- as.numeric(if (is.matrix(sl$alpha)) sl$alpha[s, ] else sl$alpha[s])
      beta_s  <- as.numeric(if (is.matrix(sl$beta))  sl$beta[s,  ] else sl$beta[s])
    } else {
      alpha_s <- as.numeric(sl$alpha[s, ])   # vector de longitud n_seg_loc
      beta_s  <- as.numeric(sl$beta[s,  ])
    }
    
    eta_s   <- sl$n[s] * sl$delta[s]
    delta_s <- sl$delta[s]
    
    tau_s <- if (n_cp == 0) {
      numeric(0)
    } else if (is.matrix(sl$tau)) {
      as.numeric(sl$tau[s, ])
    } else {
      as.numeric(sl$tau[s])
    }
    
    # Evaluar M(t_i) para todos los eventos usando la funcion con CP y
    # correccion de continuidad
    m_matrix[s, ] <- calc_M_hawkes_cp(
      t_eval    = d,
      d_hist    = d,
      alpha_vec = alpha_s,
      beta_vec  = beta_s,
      tau_vec   = tau_s,
      eta       = eta_s,
      delta     = delta_s
    )
  }
  
  # Media posterior de M(t_i) y SDM
  m_mean <- colMeans(m_matrix)
  round(mean(abs(acum_eventos - m_mean)), 4)
}

cat("Calculando SDM Hawkes (corregido con calc_M_hawkes_cp)...\n")

SDM_mat <- matrix(NA, 4, 4)
for (cp in 1:4) {
  n_cp_loc <- cp - 1
  for (mod in 1:4) {
    d  <- di[[mod]]
    sl <- samples_list[[cp]][[mod]]$BUGSoutput$sims.list
    cat(sprintf("  SDM Hawkes %dCP | mod %d...\n", n_cp_loc, mod))
    SDM_mat[mod, cp] <- calc_sdm_hawkes(d, sl, n_cp_loc)
  }
}

{
  header <- "Dataset & 0\\,CP & 1\\,CP & 2\\,CP & 3\\,CP \\\\"
  body   <- sapply(1:4, function(mod)
    d_labels_tex[mod] %+% " & " %+%
      paste(sprintf("%.4f", SDM_mat[mod, ]), collapse = " & ") %+% " \\\\")
  
  write_tabular(
    caption = "SDM por modelo y n\\'{u}mero de puntos de cambio --- Hawkes M1",
    label   = "tab:hawkes_sdm",
    col_def = "lrrrr",
    header  = header,
    body    = body,
    path    = base_dir %+% "/SDM_hawkes_todos.tex"
  )
  
  for (cp in 1:4) {
    ncp <- cp - 1
    body_cp <- sapply(1:4, function(mod)
      d_labels_tex[mod] %+% " & " %+% sprintf("%.4f", SDM_mat[mod, cp]) %+% " \\\\")
    write_tabular(
      caption = "SDM --- Hawkes M1 " %+% ncp %+% "CP",
      label   = "tab:hawkes_sdm_" %+% ncp %+% "cp",
      col_def = "lr",
      header  = "Dataset & SDM (" %+% ncp %+% "\\,CP) \\\\",
      body    = body_cp,
      path    = cp_dirs[cp] %+% "/SDM_hawkes_" %+% ncp %+% "cp.tex"
    )
  }
}

# =============================================================================
# 4. Tabla combinada WAIC + DIC + SDM
# =============================================================================

build_combined_tables <- function(WAIC_mat, DIC_mat, SDM_mat) {
  
  header_comb <- "Dataset & CP & DIC & WAIC & SDM \\\\"
  
  body_comb <- character(0)
  for (mod in 1:4) {
    first <- TRUE
    for (cp in 1:4) {
      ncp   <- cp - 1
      mc    <- if (first) d_labels_tex[mod] else ""
      first <- FALSE
      body_comb <- c(body_comb,
                     mc %+% " & " %+% ncp %+% "\\,CP & " %+%
                       sprintf("%.2f",  DIC_mat[mod,  cp]) %+% " & " %+%
                       sprintf("%.2f",  WAIC_mat[mod, cp]) %+% " & " %+%
                       sprintf("%.4f",  SDM_mat[mod,  cp]) %+% " \\\\"
      )
    }
    body_comb <- c(body_comb, "\\midrule")
  }
  body_comb <- body_comb[-length(body_comb)]
  
  write_tabular(
    caption = "Criterios de comparaci\\'{o}n de modelos: WAIC, DIC y SDM --- Hawkes M1",
    label   = "tab:hawkes_criterios",
    col_def = "llrrr",
    header  = header_comb,
    body    = body_comb,
    path    = base_dir %+% "/criterios_hawkes_todos.tex"
  )
  
  for (cp in 1:4) {
    ncp      <- cp - 1
    body_cp  <- sapply(1:4, function(mod)
      d_labels_tex[mod] %+% " & " %+%
        sprintf("%.2f",  DIC_mat[mod,  cp]) %+% " & " %+%
        sprintf("%.2f",  WAIC_mat[mod, cp]) %+% " & " %+%
        sprintf("%.4f",  SDM_mat[mod,  cp]) %+% " \\\\"
    )
    write_tabular(
      caption = "Criterios de comparaci\\'{o}n: WAIC, DIC y SDM --- Hawkes M1 " %+%
        ncp %+% "CP",
      label   = "tab:hawkes_criterios_" %+% ncp %+% "cp",
      col_def = "lrrr",
      header  = "Dataset & DIC & WAIC & SDM \\\\",
      body    = body_cp,
      path    = cp_dirs[cp] %+% "/criterios_hawkes_" %+% ncp %+% "cp.tex"
    )
  }
}

build_combined_tables(WAIC_mat, DIC_mat, SDM_mat)

# =============================================================================
# 6. Tabla vida media del kernel: log(2)/delta
# =============================================================================

VM_mat <- matrix(NA, 4, 4)
for (cp in 1:4)
  for (mod in 1:4) {
    su  <- samples_list[[cp]][[mod]]$BUGSoutput$summary
    dlt <- su[grep("^delta$", rownames(su)), "mean"]
    VM_mat[mod, cp] <- round(log(2) / dlt, 2)
  }

{
  header <- "Dataset & 0\\,CP & 1\\,CP & 2\\,CP & 3\\,CP \\\\"
  body   <- sapply(1:4, function(mod)
    d_labels_tex[mod] %+% " & " %+%
      paste(sprintf("%.2f", VM_mat[mod, ]), collapse = " & ") %+% " \\\\")
  
  write_tabular(
    caption  = "Vida media del kernel $\\log(2)/\\delta$ (d\\'{i}as) --- Hawkes M1",
    label    = "tab:hawkes_halflife",
    col_def  = "lrrrr",
    header   = header,
    body     = body,
    footnote = "La vida media mide cu\\'{a}nto tarda el efecto de auto-excitaci\\'{o}n en reducirse a la mitad.",
    path     = base_dir %+% "/vida_media_hawkes_todos.tex"
  )
  
  for (cp in 1:4) {
    ncp <- cp - 1
    body_cp <- sapply(1:4, function(mod)
      d_labels_tex[mod] %+% " & " %+% sprintf("%.2f", VM_mat[mod, cp]) %+% " \\\\")
    write_tabular(
      caption = "Vida media del kernel --- Hawkes M1 " %+% ncp %+% "CP",
      label   = "tab:hawkes_halflife_" %+% ncp %+% "cp",
      col_def = "lr",
      header  = "Dataset & $\\log(2)/\\delta$ (d\\'{i}as) \\\\",
      body    = body_cp,
      path    = cp_dirs[cp] %+% "/vida_media_hawkes_" %+% ncp %+% "cp.tex"
    )
  }
}