# =============================================================================
# Tablas de resultados — NHPP Weibull
# Requiere: samples_NHPP0, samples_NHPP1, samples_NHPP2, samples_NHPP3
# Produce tablas .tex con convencion booktabs/longtable del articulo base:
#   1. Tabla de parametros (Model | Parameters | Prior | Mean | SD | CI 95%)
#   2. Tabla DIC comparativa (0CP … 3CP)
#   3. Tabla WAIC comparativa
#   4. Tabla SDM comparativa
# =============================================================================

library(loo)

load("NHPP_0cp.RData")
load("NHPP_1cp.RData")
load("NHPP_2cp.RData")
load("NHPP_3cp.RData")

samples_NHPP0 <- samples_resultados0

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
samples_list <- list(samples_NHPP0, samples_NHPP1,
                     samples_NHPP2, samples_NHPP3)

prior_alpha <- "$\\text{Unif}(0,\\,3)$"
prior_beta  <- "$\\text{Unif}(0,\\,100)$"
prior_tau   <- c("$\\text{Unif}(0,\\,365)$",
                 "$\\text{Unif}(365,\\,730)$",
                 "$\\text{Unif}(730,\\,1096)$")

# Carpetas
base_dir <- "10_resultados_tabulares"
cp_dirs  <- base_dir %+% "/resultados_" %+% 0:3 %+% "cp"
for (d in c(base_dir, cp_dirs)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Cuenta columnas en col_def (solo letras l, r, c)
ncols_of <- function(col_def) nchar(gsub("[^lrc]", "", col_def))

# Escribe tabla tabular simple (table flotante)
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

# Escribe tabla longtable (para tablas largas de parametros)
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
#    Columnas: Model | Parameters | Prior Distribution | Mean | SD | CI 95%
# =============================================================================

n_seg <- c(1, 2, 3, 4)

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
      # Fallback: si no tiene index (modelo 0CP)
      if (nrow(a_rows) == 0) a_rows <- su[grep("^alpha$", rownames(su)), , drop = FALSE]
      if (nrow(b_rows) == 0) b_rows <- su[grep("^beta$",  rownames(su)), , drop = FALSE]
      
      mc <- if (first) nombre else ""
      body <- c(body, emit(mc, "$\\alpha_{" %+% i %+% "}$", prior_alpha, a_rows))
      first <- FALSE
      body <- c(body, emit("",  "$\\beta_{" %+% i %+% "}$",  prior_beta,  b_rows))
    }
    
    if (ncp > 0) {
      for (j in 1:ncp) {
        t_rows <- su[grep("^tau\\[" %+% j %+% "\\]", rownames(su)), , drop = FALSE]
        if (nrow(t_rows) == 0) t_rows <- su[grep("^tau$", rownames(su)), , drop = FALSE]
        body <- c(body, emit("", "$\\tau_{" %+% j %+% "}$", prior_tau[j], t_rows))
      }
    }
    
    body <- c(body, "\\midrule")
  }
  body <- body[-length(body)]  # quitar ultimo \midrule
  
  col_def <- "llllrrr"
  header  <- paste("Model & Parameters & Prior Distribution &",
                   "Mean & SD & Credible intervals 95\\% \\\\")
  caption <- "Prior distributions and estimated parameters --- NHPP Weibull " %+% ncp %+% "CP"
  label   <- "tab:nhpp_params_" %+% ncp %+% "cp"
  path    <- cp_dirs[cp] %+% "/parametros_NHPP_" %+% ncp %+% "cp.tex"
  
  write_longtable(caption, label, col_def, header, body, path)
}

# =============================================================================
# 2. Tabla DIC
# =============================================================================

DIC_mat <- matrix(NA, 4, 4)
for (cp in 1:4)
  for (mod in 1:4)
    DIC_mat[mod, cp] <- round(samples_list[[cp]][[mod]]$BUGSoutput$DIC, 0)

{
  header <- "Dataset & 0\\,CP & 1\\,CP & 2\\,CP & 3\\,CP \\\\"
  body   <- sapply(1:4, function(mod)
    d_labels_tex[mod] %+% " & " %+%
      paste(DIC_mat[mod, ], collapse = " & ") %+% " \\\\")
  
  write_tabular(
    caption  = "DIC por modelo y n\\'{u}mero de puntos de cambio --- NHPP Weibull",
    label    = "tab:nhpp_dic",
    col_def  = "lrrrr",
    header   = header,
    body     = body,
    path     = base_dir %+% "/DIC_NHPP_todos.tex"
  )
  
  for (cp in 1:4) {
    ncp  <- cp - 1
    body_cp <- sapply(1:4, function(mod)
      d_labels_tex[mod] %+% " & " %+% DIC_mat[mod, cp] %+% " \\\\")
    write_tabular(
      caption = "DIC --- NHPP Weibull " %+% ncp %+% "CP",
      label   = "tab:nhpp_dic_" %+% ncp %+% "cp",
      col_def = "lr",
      header  = "Dataset & DIC (" %+% ncp %+% "\\,CP) \\\\",
      body    = body_cp,
      path    = cp_dirs[cp] %+% "/DIC_NHPP_" %+% ncp %+% "cp.tex"
    )
  }
    }

# =============================================================================
# 3. Tabla WAIC
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
    caption = "WAIC por modelo y n\\'{u}mero de puntos de cambio --- NHPP Weibull",
    label   = "tab:nhpp_waic",
    col_def = "lrrrr",
    header  = header,
    body    = body,
    path    = base_dir %+% "/WAIC_NHPP_todos.tex"
  )
  
  for (cp in 1:4) {
    ncp <- cp - 1
    body_cp <- sapply(1:4, function(mod)
      d_labels_tex[mod] %+% " & " %+% sprintf("%.2f", WAIC_mat[mod, cp]) %+% " \\\\")
    write_tabular(
      caption = "WAIC --- NHPP Weibull " %+% ncp %+% "CP",
      label   = "tab:nhpp_waic_" %+% ncp %+% "cp",
      col_def = "lr",
      header  = "Dataset & WAIC (" %+% ncp %+% "\\,CP) \\\\",
      body    = body_cp,
      path    = cp_dirs[cp] %+% "/WAIC_NHPP_" %+% ncp %+% "cp.tex"
    )
  }
}

# =============================================================================
# 4. Tabla SDM
# =============================================================================

SDM_mat <- matrix(NA, 4, 4)
for (cp in 1:4)
  for (mod in 1:4) {
    d    <- di[[mod]]
    acum <- sapply(1:t_total, function(i) sum(d <= i))
    su   <- samples_list[[cp]][[mod]]$BUGSoutput$summary
    mm   <- su[grep("^m", rownames(su)), "mean"]
    SDM_mat[mod, cp] <- round(sum(abs(acum[d] - mm)) / length(d), 4)
  }

{
  header <- "Dataset & 0\\,CP & 1\\,CP & 2\\,CP & 3\\,CP \\\\"
  body   <- sapply(1:4, function(mod)
    d_labels_tex[mod] %+% " & " %+%
      paste(sprintf("%.4f", SDM_mat[mod, ]), collapse = " & ") %+% " \\\\")
  
  write_tabular(
    caption = "SDM por modelo y n\\'{u}mero de puntos de cambio --- NHPP Weibull",
    label   = "tab:nhpp_sdm",
    col_def = "lrrrr",
    header  = header,
    body    = body,
    path    = base_dir %+% "/SDM_NHPP_todos.tex"
  )
  
  for (cp in 1:4) {
    ncp <- cp - 1
    body_cp <- sapply(1:4, function(mod)
      d_labels_tex[mod] %+% " & " %+% sprintf("%.4f", SDM_mat[mod, cp]) %+% " \\\\")
    write_tabular(
      caption = "SDM --- NHPP Weibull " %+% ncp %+% "CP",
      label   = "tab:nhpp_sdm_" %+% ncp %+% "cp",
      col_def = "lr",
      header  = "Dataset & SDM (" %+% ncp %+% "\\,CP) \\\\",
      body    = body_cp,
      path    = cp_dirs[cp] %+% "/SDM_NHPP_" %+% ncp %+% "cp.tex"
    )
  }
}