# =============================================================================
# Tablas de resultados — Hawkes M1
# Requiere: samples_hawkes0, samples_hawkes1, samples_hawkes2, samples_hawkes3
# Produce tablas .tex con convencion booktabs/longtable del articulo base:
#   1. Tabla de parametros (Model | Parameters | Prior | Mean | SD | CI 95%)
#      alpha/beta por segmento + tau + n, delta, eta globales
#   2. Tabla WAIC (DIC excluido — incompatible con NHPP)
#   3. Tabla SDM
#   4. Tabla vida media del kernel log(2)/delta
# =============================================================================

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
    
    # Parametros globales del kernel Hawkes
    for (par_info in list(
      list(pat = "^n$",     tex = "$n$",      prior = prior_n),
      list(pat = "^delta$", tex = "$\\delta$", prior = prior_delta),
      list(pat = "^eta$",   tex = "$\\eta$",   prior = "derivado: $n \\cdot \\delta$"))) {
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
# 2. Tabla WAIC (DIC excluido — incompatible por zeros trick)
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
    caption  = "WAIC por modelo y n\\'{u}mero de puntos de cambio --- Hawkes M1",
    label    = "tab:hawkes_waic",
    col_def  = "lrrrr",
    header   = header,
    body     = body,
    footnote = "El DIC no se reporta para el proceso de Hawkes por incompatibilidad con el NHPP debida a la constante $C$ del \\textit{zeros trick} (diferencia sistem\\'{a}tica de $2KC$ en la devianza).",
    path     = base_dir %+% "/WAIC_hawkes_todos.tex"
  )
  
  for (cp in 1:4) {
    ncp <- cp - 1
    body_cp <- sapply(1:4, function(mod)
      d_labels_tex[mod] %+% " & " %+% sprintf("%.2f", WAIC_mat[mod, cp]) %+% " \\\\")
    write_tabular(
      caption  = "WAIC --- Hawkes M1 " %+% ncp %+% "CP",
      label    = "tab:hawkes_waic_" %+% ncp %+% "cp",
      col_def  = "lr",
      header   = "Dataset & WAIC (" %+% ncp %+% "\\,CP) \\\\",
      body     = body_cp,
      path     = cp_dirs[cp] %+% "/WAIC_hawkes_" %+% ncp %+% "cp.tex"
    )
  }
}

# =============================================================================
# 3. Tabla SDM — calculo externo via sims.list (M(t) forma cerrada)
# M(t) = (t/beta)^alpha + (eta/delta) * sum_{d_k<t}(1 - exp(-delta*(t-d_k)))
# Se propaga la incertidumbre MCMC y se toma la media posterior de M(t_i)
# =============================================================================

# Funcion auxiliar: media posterior de M(t) en los tiempos de evento
# Itera sobre n_sims muestras MCMC para propagar incertidumbre
calc_sdm_hawkes <- function(d, sl, n_seg_mod) {
  # sl: sims.list del modelo
  # n_seg_mod: numero de segmentos (1 para 0CP, 2 para 1CP, etc.)
  n_sims <- length(sl$n)
  acum   <- sapply(1:t_total, function(i) sum(d <= i))
  
  # Matriz: filas = simulaciones, columnas = tiempos de evento d
  m_matrix <- matrix(NA, nrow = n_sims, ncol = length(d))
  
  for (s in 1:n_sims) {
    # Extraer alpha y beta de la simulacion s
    # Para 0CP son vectores; para CP>0 son matrices [sim, segmento]
    if (n_seg_mod == 1) {
      alpha_s <- if (is.matrix(sl$alpha)) sl$alpha[s, 1] else sl$alpha[s]
      beta_s  <- if (is.matrix(sl$beta))  sl$beta[s,  1] else sl$beta[s]
    } else {
      # Para modelos con CP usamos el alpha/beta del ultimo segmento
      # como aproximacion global para el calculo del SDM externo.
      # El compensador exacto por regimenes requeriria integrar por tramos,
      # lo que es computacionalmente costoso. Esta aproximacion es coherente
      # con el uso de la media posterior de alpha y beta en el script de graficas.
      alpha_s <- sl$alpha[s, n_seg_mod]
      beta_s  <- sl$beta[s,  n_seg_mod]
    }
    eta_s   <- sl$n[s] * sl$delta[s]
    delta_s <- sl$delta[s]
    
    m_matrix[s, ] <- sapply(d, function(tt) {
      prev <- d[d < tt]
      (tt / beta_s)^alpha_s +
        if (length(prev) > 0)
          (eta_s / delta_s) * sum(1 - exp(-delta_s * (tt - prev)))
      else 0
    })
  }
  
  # Media posterior de M(t_i) para cada tiempo de evento
  mm_mean <- colMeans(m_matrix)
  
  # SDM = (1/K) * sum |N(t_i) - M_hat(t_i)|
  round(sum(abs(acum[d] - mm_mean)) / length(d), 4)
}

cat("Calculando SDM Hawkes (metodo externo via sims.list)...\n")
cat("Esto puede tardar varios minutos por modelo.\n")

SDM_mat <- matrix(NA, 4, 4)
for (cp in 1:4) {
  n_seg_mod <- cp  # 0CP -> 1 seg, 1CP -> 2 seg, etc.
  for (mod in 1:4) {
    d  <- di[[mod]]
    sl <- samples_list[[cp]][[mod]]$BUGSoutput$sims.list
    cat(sprintf("  SDM Hawkes %dCP | mod %d...\n", cp - 1, mod))
    SDM_mat[mod, cp] <- calc_sdm_hawkes(d, sl, n_seg_mod)
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
# 4. Tabla vida media del kernel: log(2)/delta
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