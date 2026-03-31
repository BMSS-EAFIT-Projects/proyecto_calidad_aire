# =============================================================================
# Tablas comparativas — NHPP vs Hawkes M1
# Requiere: todos los RData de NHPP y Hawkes cargados
# Metricas compatibles: WAIC y SDM
# DIC EXCLUIDO — incompatible por diferencia sistematica de 2KC
#   en la devianza (constante C del zeros trick de Hawkes)
# Produce tablas .tex con convencion booktabs del articulo base
# =============================================================================

library(loo)

load("NHPP_0cp.RData");   load("NHPP_1cp.RData")
load("NHPP_2cp.RData");   load("NHPP_3cp.RData")
load("hawkes_0cp.RData"); load("hawkes_1cp.RData")
load("hawkes_2cp.RData"); load("hawkes_3cp.RData")

`%+%` <- paste0

d_labels_tex <- c(
  "PM$_{2.5}$ Colombia 37 --- Bogot\\'{a}",
  "PM$_{2.5}$ Colombia 37 --- Medell\\'{i}n",
  "PM$_{2.5}$ WHO 25 --- Bogot\\'{a}",
  "PM$_{2.5}$ WHO 25 --- Medell\\'{i}n"
)

di <- list(excePM25col_Bog, excePM25col_Med,
           excePM25who_Bog, excePM25who_Med)

t_total     <- 1096
nhpp_list   <- list(samples_NHPP0,   samples_NHPP1,   samples_NHPP2,   samples_NHPP3)
hawkes_list <- list(samples_hawkes0, samples_hawkes1, samples_hawkes2, samples_hawkes3)
cp_labels   <- c("0\\,CP", "1\\,CP", "2\\,CP", "3\\,CP")

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

# =============================================================================
# 1. WAIC comparativo NHPP vs Hawkes — una tabla por CP + tabla completa
# =============================================================================

WAIC_nhpp   <- matrix(NA, 4, 4)
WAIC_hawkes <- matrix(NA, 4, 4)

for (cp in 1:4)
  for (mod in 1:4) {
    ll_n <- nhpp_list[[cp]][[mod]]$BUGSoutput$sims.list$phi2
    ll_h <- hawkes_list[[cp]][[mod]]$BUGSoutput$sims.list$phi2
    WAIC_nhpp[mod, cp]   <- round(waic(ll_n)$estimates["elpd_waic","Estimate"] * -2, 2)
    WAIC_hawkes[mod, cp] <- round(waic(ll_h)$estimates["elpd_waic","Estimate"] * -2, 2)
  }

# Tabla completa (8 columnas: NHPP_0CP … NHPP_3CP | Hawkes_0CP … Hawkes_3CP)
{
  header <- paste(
    "Dataset",
    paste(sapply(cp_labels, function(l) "NHPP " %+% l), collapse = " & "),
    paste(sapply(cp_labels, function(l) "Hawkes " %+% l), collapse = " & "),
    sep = " & "
  ) %+% " \\\\"
  
  body <- sapply(1:4, function(mod)
    d_labels_tex[mod] %+% " & " %+%
      paste(sprintf("%.2f", WAIC_nhpp[mod, ]),   collapse = " & ") %+% " & " %+%
      paste(sprintf("%.2f", WAIC_hawkes[mod, ]), collapse = " & ") %+% " \\\\")
  
  write_tabular(
    caption  = "WAIC comparativo NHPP vs Hawkes M1 por n\\'{u}mero de puntos de cambio",
    label    = "tab:comp_waic",
    col_def  = "lrrrrrrrr",
    header   = header,
    body     = body,
    footnote = "El DIC no se reporta para el modelo Hawkes por incompatibilidad debida al \\textit{zeros trick}.",
    path     = base_dir %+% "/WAIC_comparativo_NHPP_vs_Hawkes.tex"
  )
}

# Tabla por CP (columnas: NHPP | Hawkes para ese numero de CP)
for (cp in 1:4) {
  ncp    <- cp - 1
  header <- "Dataset & NHPP " %+% cp_labels[cp] %+% " & Hawkes " %+% cp_labels[cp] %+% " \\\\"
  body   <- sapply(1:4, function(mod)
    d_labels_tex[mod] %+% " & " %+%
      sprintf("%.2f", WAIC_nhpp[mod, cp])   %+% " & " %+%
      sprintf("%.2f", WAIC_hawkes[mod, cp]) %+% " \\\\")
  
  write_tabular(
    caption = "WAIC comparativo --- NHPP vs Hawkes M1, " %+% ncp %+% " punto(s) de cambio",
    label   = "tab:comp_waic_" %+% ncp %+% "cp",
    col_def = "lrr",
    header  = header,
    body    = body,
    path    = cp_dirs[cp] %+% "/WAIC_comparativo_" %+% ncp %+% "cp.tex"
  )
}

# =============================================================================
# 2. SDM comparativo NHPP vs Hawkes
# SDM NHPP   : desde summary de JAGS (m[i] monitoreado directamente)
# SDM Hawkes : calculo externo via sims.list (M(t) forma cerrada)
#   M(t) = (t/beta)^alpha + (eta/delta)*sum_{d_k<t}(1-exp(-delta*(t-d_k)))
# =============================================================================

# Funcion auxiliar: SDM Hawkes externo propagando muestras MCMC
calc_sdm_hawkes_comp <- function(d, sl, n_seg_mod) {
  n_sims <- length(sl$n)
  acum   <- sapply(1:t_total, function(i) sum(d <= i))
  
  m_matrix <- matrix(NA, nrow = n_sims, ncol = length(d))
  
  for (s in 1:n_sims) {
    if (n_seg_mod == 1) {
      alpha_s <- if (is.matrix(sl$alpha)) sl$alpha[s, 1] else sl$alpha[s]
      beta_s  <- if (is.matrix(sl$beta))  sl$beta[s,  1] else sl$beta[s]
    } else {
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
  
  mm_mean <- colMeans(m_matrix)
  round(sum(abs(acum[d] - mm_mean)) / length(d), 4)
}

SDM_nhpp   <- matrix(NA, 4, 4)
SDM_hawkes <- matrix(NA, 4, 4)

cat("Calculando SDM comparativo (Hawkes via sims.list)...\n")

for (cp in 1:4)
  for (mod in 1:4) {
    d    <- di[[mod]]
    acum <- sapply(1:t_total, function(i) sum(d <= i))
    
    # SDM NHPP — desde summary de JAGS
    su_n  <- nhpp_list[[cp]][[mod]]$BUGSoutput$summary
    mm_n  <- su_n[grep("^m", rownames(su_n)), "mean"]
    SDM_nhpp[mod, cp] <- round(sum(abs(acum[d] - mm_n)) / length(d), 4)
    
    # SDM Hawkes — calculo externo
    cat(sprintf("  SDM Hawkes %dCP | mod %d...\n", cp - 1, mod))
    sl_h <- hawkes_list[[cp]][[mod]]$BUGSoutput$sims.list
    SDM_hawkes[mod, cp] <- calc_sdm_hawkes_comp(d, sl_h, n_seg_mod = cp)
  }

# Tabla completa
{
  header <- paste(
    "Dataset",
    paste(sapply(cp_labels, function(l) "NHPP " %+% l), collapse = " & "),
    paste(sapply(cp_labels, function(l) "Hawkes " %+% l), collapse = " & "),
    sep = " & "
  ) %+% " \\\\"
  
  fmt4v <- function(v) sapply(v, function(x) if (is.na(x)) "---" else sprintf("%.4f", x))
  
  body <- sapply(1:4, function(mod)
    d_labels_tex[mod] %+% " & " %+%
      paste(fmt4v(SDM_nhpp[mod, ]),   collapse = " & ") %+% " & " %+%
      paste(fmt4v(SDM_hawkes[mod, ]), collapse = " & ") %+% " \\\\")
  
  write_tabular(
    caption = "SDM comparativo NHPP vs Hawkes M1 por n\\'{u}mero de puntos de cambio",
    label   = "tab:comp_sdm",
    col_def = "lrrrrrrrr",
    header  = header,
    body    = body,
    path    = base_dir %+% "/SDM_comparativo_NHPP_vs_Hawkes.tex"
  )
}

# Tabla por CP
fmt4s <- function(x) if (is.na(x)) "---" else sprintf("%.4f", x)

for (cp in 1:4) {
  ncp    <- cp - 1
  header <- "Dataset & NHPP " %+% cp_labels[cp] %+% " & Hawkes " %+% cp_labels[cp] %+% " \\\\"
  body   <- sapply(1:4, function(mod)
    d_labels_tex[mod] %+% " & " %+%
      fmt4s(SDM_nhpp[mod, cp])   %+% " & " %+%
      fmt4s(SDM_hawkes[mod, cp]) %+% " \\\\")
  
  write_tabular(
    caption = "SDM comparativo --- NHPP vs Hawkes M1, " %+% ncp %+% " punto(s) de cambio",
    label   = "tab:comp_sdm_" %+% ncp %+% "cp",
    col_def = "lrr",
    header  = header,
    body    = body,
    path    = cp_dirs[cp] %+% "/SDM_comparativo_" %+% ncp %+% "cp.tex"
  )
}

# =============================================================================
# 3. Tabla resumen: mejor modelo por contaminante segun WAIC y SDM
# =============================================================================

{
  best_nhpp_waic   <- apply(WAIC_nhpp,   1, which.min)
  best_hawkes_waic <- apply(WAIC_hawkes, 1, which.min)
  best_nhpp_sdm    <- apply(SDM_nhpp,    1, which.min)
  best_hawkes_sdm  <- apply(SDM_hawkes,  1, function(x) which.min(replace(x, is.na(x), Inf)))
  
  cp_short <- c("0CP", "1CP", "2CP", "3CP")
  
  header <- paste("Dataset",
                  "Best NHPP", "WAIC",
                  "Best Hawkes", "WAIC",
                  "Winner (WAIC)",
                  "Best NHPP", "SDM",
                  "Best Hawkes", "SDM",
                  "Winner (SDM)",
                  sep = " & ") %+% " \\\\"
  
  # Formatea un valor numerico: NA -> "---"
  fmt2 <- function(x) if (is.na(x)) "---" else sprintf("%.2f", x)
  fmt4 <- function(x) if (is.na(x)) "---" else sprintf("%.4f", x)
  
  body <- sapply(1:4, function(mod) {
    bn_w <- best_nhpp_waic[mod];   bh_w <- best_hawkes_waic[mod]
    bn_s <- best_nhpp_sdm[mod];    bh_s <- best_hawkes_sdm[mod]
    
    waic_n <- WAIC_nhpp[mod, bn_w];   waic_h <- WAIC_hawkes[mod, bh_w]
    sdm_n  <- SDM_nhpp[mod,  bn_s];   sdm_h  <- SDM_hawkes[mod,  bh_s]
    
    # isTRUE protege contra NA en la comparacion
    win_w <- if (isTRUE(waic_n < waic_h)) "NHPP" else if (isTRUE(waic_h < waic_n)) "Hawkes" else "---"
    win_s <- if (isTRUE(sdm_n  < sdm_h))  "NHPP" else if (isTRUE(sdm_h  < sdm_n))  "Hawkes" else "---"
    
    paste(d_labels_tex[mod],
          cp_short[bn_w], fmt2(waic_n),
          cp_short[bh_w], fmt2(waic_h),
          "\\textbf{" %+% win_w %+% "}",
          cp_short[bn_s], fmt4(sdm_n),
          cp_short[bh_s], fmt4(sdm_h),
          "\\textbf{" %+% win_s %+% "}",
          sep = " & ") %+% " \\\\"
  })
  
  write_tabular(
    caption  = "Resumen de selecci\\'{o}n de modelos: mejor configuraci\\'{o}n por contaminante --- NHPP vs Hawkes M1",
    label    = "tab:comp_resumen",
    col_def  = "lrrrrlrrrrl",
    header   = header,
    body     = body,
    footnote = "Se reporta el n\\'{u}mero de puntos de cambio del mejor modelo dentro de cada familia y la m\\'{e}trica correspondiente. El ganador se determina comparando los mejores modelos de cada familia.",
    path     = base_dir %+% "/resumen_seleccion_NHPP_vs_Hawkes.tex"
  )
}