# =============================================================================
# Tests de Convergencia — 0 Change Points
# Corre sobre: samples_NHPP0 y samples_hawkes0
# Guarda en:   09_resultados_convergencia/resultados_convergencia_0cp/
# =============================================================================

library(xtable)
library(rjags)
library(lattice)
library(stats)
library(bayesplot)
library(loo)
library(R2jags)
library(coda)

source("07_utilidades/crear_mcmc_list.R")

load("NHPP_0cp.RData")
load("hawkes_0cp.RData")

dir_out <- "09_resultados_convergencia/resultados_convergencia_0cp"

d_names <- c("PM25col_Bog", "PM25col_Med", "PM25who_Bog", "PM25who_Med")

di <- list(excePM25col_Bog, excePM25col_Med,
           excePM25who_Bog, excePM25who_Med)

# -----------------------------------------------------------------------------
# Función auxiliar: extrae sims como mcmc.list con nombres
# tipo = "NHPP" o "Hawkes"
# -----------------------------------------------------------------------------
obtener_sims <- function(modelo, samples_list, tipo) {
  sl <- samples_list[[modelo]]$BUGSoutput$sims.list
  nc <- samples_list[[modelo]]$BUGSoutput$n.chains

  alphas <- crear_mcmc_list(as.matrix(sl$alpha), nc)
  betas  <- crear_mcmc_list(as.matrix(sl$beta),  nc)

  out <- list(alphas = alphas, betas = betas)

  if (tipo == "Hawkes") {
    out$ns     <- crear_mcmc_list(as.matrix(sl$n),     nc)
    out$deltas <- crear_mcmc_list(as.matrix(sl$delta), nc)
  }

  return(out)
}

# -----------------------------------------------------------------------------
# Función auxiliar: tests numéricos para una cadena mcmc.list
# Devuelve lista con raftery I, heidel p-valores, gelman upper CI, geweke z
# -----------------------------------------------------------------------------
procesar_cadena <- function(cadena) {
  cadena_combinada <- mcmc(do.call(rbind, cadena))
  list(
    raftery = raftery.diag(cadena_combinada)$resmatrix[, "I"],
    heidel1 = heidel.diag(cadena)[[1]][, 3],
    heidel2 = heidel.diag(cadena)[[2]][, 3],
    gelman  = gelman.diag(cadena)[1]$psrf[, 2],
    geweke1 = geweke.diag(cadena)[[1]]$z,
    geweke2 = geweke.diag(cadena)[[2]]$z
  )
}

# -----------------------------------------------------------------------------
# Función auxiliar: genera filas LaTeX para la tabla de convergencia
# -----------------------------------------------------------------------------
generar_filas <- function(res, nombre_param, n_params) {
  filas <- ""
  for (i in 1:n_params) {
    fila <- sprintf(
      "\\multirow{2}{*}{$\\%s_%d$} & Chain 1 & \\multirow{2}{*}{%.2f} & %.3f & \\multirow{2}{*}{%.3f} & %.2f \\\\\n & Chain 2 & & %.3f & & %.2f \\\\\n",
      nombre_param, i,
      res$raftery[i],
      res$heidel1[i], res$gelman[i], res$geweke1[i],
      res$heidel2[i], res$geweke2[i]
    )
    filas <- paste0(filas, fila)
  }
  return(filas)
}

# =============================================================================
# LOOP PRINCIPAL — corre sobre los 4 contaminantes x 2 modelos
# =============================================================================

for (mod in 1:4) {

  nombre <- d_names[mod]
  cat(sprintf("\n====== Convergencia | 0CP | %s ======\n", nombre))

  for (tipo in c("NHPP", "Hawkes")) {

    samples_list <- if (tipo == "NHPP") samples_resultados0 else samples_hawkes0
    sims <- obtener_sims(mod, samples_list, tipo)

    prefijo <- sprintf("%s_%s_0cp", tipo, nombre)

    # -------------------------------------------------------------------------
    # 1. Density plots
    # -------------------------------------------------------------------------
    png(file.path(dir_out, paste0("density_alpha_", prefijo, ".png")),
        width = 1800, height = 300, res = 120)
    print(densityplot(sims$alphas, ylab = "", main = paste("Density alpha |", tipo, "|", nombre),
                      col = c("red", "#2CA02C"), layout = c(1, 1)))
    dev.off()

    png(file.path(dir_out, paste0("density_beta_", prefijo, ".png")),
        width = 1800, height = 300, res = 120)
    print(densityplot(sims$betas, ylab = "", main = paste("Density beta |", tipo, "|", nombre),
                      col = c("red", "#2CA02C"), layout = c(1, 1)))
    dev.off()

    if (tipo == "Hawkes") {
      png(file.path(dir_out, paste0("density_n_", prefijo, ".png")),
          width = 1800, height = 300, res = 120)
      print(densityplot(sims$ns, ylab = "", main = paste("Density n |", nombre),
                        col = c("red", "#2CA02C"), layout = c(1, 1)))
      dev.off()

      png(file.path(dir_out, paste0("density_delta_", prefijo, ".png")),
          width = 1800, height = 300, res = 120)
      print(densityplot(sims$deltas, ylab = "", main = paste("Density delta |", nombre),
                        col = c("red", "#2CA02C"), layout = c(1, 1)))
      dev.off()
    }

    # -------------------------------------------------------------------------
    # 2. Trace plots
    # -------------------------------------------------------------------------
    png(file.path(dir_out, paste0("trace_alpha_", prefijo, ".png")),
        width = 1800, height = 300, res = 120)
    print(xyplot(sims$alphas, main = paste("Trace alpha |", tipo, "|", nombre),
                 col = c("red", "#2CA02C"), layout = c(1, 1)))
    dev.off()

    png(file.path(dir_out, paste0("trace_beta_", prefijo, ".png")),
        width = 1800, height = 300, res = 120)
    print(xyplot(sims$betas, main = paste("Trace beta |", tipo, "|", nombre),
                 col = c("red", "#2CA02C"), layout = c(1, 1)))
    dev.off()

    if (tipo == "Hawkes") {
      png(file.path(dir_out, paste0("trace_n_", prefijo, ".png")),
          width = 1800, height = 300, res = 120)
      print(xyplot(sims$ns, main = paste("Trace n |", nombre),
                   col = c("red", "#2CA02C"), layout = c(1, 1)))
      dev.off()

      png(file.path(dir_out, paste0("trace_delta_", prefijo, ".png")),
          width = 1800, height = 300, res = 120)
      print(xyplot(sims$deltas, main = paste("Trace delta |", nombre),
                   col = c("red", "#2CA02C"), layout = c(1, 1)))
      dev.off()
    }

    # -------------------------------------------------------------------------
    # 3. Autocorrelation plots
    # -------------------------------------------------------------------------
    png(file.path(dir_out, paste0("autocorr_alpha_", prefijo, ".png")),
        width = 1800, height = 300, res = 120)
    par(mfcol = c(1, 2))
    autocorr.plot(sims$alphas, auto.layout = FALSE)
    dev.off()

    png(file.path(dir_out, paste0("autocorr_beta_", prefijo, ".png")),
        width = 1800, height = 300, res = 120)
    par(mfcol = c(1, 2))
    autocorr.plot(sims$betas, auto.layout = FALSE)
    dev.off()

    if (tipo == "Hawkes") {
      png(file.path(dir_out, paste0("autocorr_n_", prefijo, ".png")),
          width = 1800, height = 300, res = 120)
      par(mfcol = c(1, 2))
      autocorr.plot(sims$ns, auto.layout = FALSE)
      dev.off()

      png(file.path(dir_out, paste0("autocorr_delta_", prefijo, ".png")),
          width = 1800, height = 300, res = 120)
      par(mfcol = c(1, 2))
      autocorr.plot(sims$deltas, auto.layout = FALSE)
      dev.off()
    }

    # -------------------------------------------------------------------------
    # 4. Tests numericos + tabla LaTeX
    # -------------------------------------------------------------------------
    res_alpha <- procesar_cadena(sims$alphas)
    res_beta  <- procesar_cadena(sims$betas)

    filas <- paste0(
      generar_filas(res_alpha, "alpha", 1),
      "\\midrule\n",
      generar_filas(res_beta,  "beta",  1)
    )

    if (tipo == "Hawkes") {
      res_n     <- procesar_cadena(sims$ns)
      res_delta <- procesar_cadena(sims$deltas)
      filas <- paste0(filas,
        "\\midrule\n",
        generar_filas(res_n,     "n",     1),
        "\\midrule\n",
        generar_filas(res_delta, "delta", 1)
      )
    }

    tabla_latex <- paste0(
      "\\begin{table}[H]\n\\centering\n",
      sprintf("\\caption{Convergence diagnostics — %s 0CP — %s}\n", tipo, nombre),
      "\\begin{tabular}{lccccc}\n\\toprule\n",
      "\\multirow{2}{*}{Parameter} & \\multirow{2}{*}{Chain} & Raftery-Lewis & Heidel-Welch & \\multirow{2}{*}{Gelman-Rubin} & \\multirow{2}{*}{Geweke} \\\\\n",
      " & & (I) & (p-value) & & (z) \\\\\n\\midrule\n",
      filas,
      "\\bottomrule\n\\end{tabular}\n\\end{table}"
    )

    writeLines(tabla_latex,
               file.path(dir_out, paste0("tabla_convergencia_", prefijo, ".tex")))

    cat(sprintf("  [%s] Listo — graficas y tabla guardadas\n", tipo))
  }
}

cat(sprintf("\nTodo guardado en: %s\n", dir_out))
