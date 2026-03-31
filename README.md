# Análisis Bayesiano de Calidad del Aire — Bogotá y Medellín

Modelos de puntos de cambio (Non-Homogeneous Poisson Process) para excedencias
de contaminantes atmosféricos (O₃, PM10, PM2.5) en el período 2018–2020,
estimados con JAGS a través de `R2jags`.

---

## Estructura del proyecto

```
proyecto_calidad_aire/
│
├── 01_datos/
│   ├── limpieza_Bogota.R        # Lectura y procesamiento datos RMCAB (Excel)
│   └── limpieza_Medellin.R      # Lectura y procesamiento datos Medellín (CSV)
│
├── 02_modelos_JAGS/
│   └── bug_files/
│       ├── 0_puntos_cambio.bug  # Modelo Weibull sin puntos de cambio
│       ├── 1_punto_cambio.bug   # Modelo con 1 punto de cambio
│       ├── 2_puntos_cambio.bug  # Modelo con 2 puntos de cambio
│       └── 3_puntos_cambio.bug  # Modelo con 3 puntos de cambio
│
├── 03_estimacion/
│   ├── runJAGS_0cp_V2.R         # Estimación 0 CP — todos los modelos (versión final)
│   ├── runJAGS_1cp_V2.R         # Estimación 1 CP — todos los modelos (versión final)
│   ├── runJAGS_2cp_V2.R         # Estimación 2 CP — todos los modelos (versión final)
│   ├── runJAGS_3cp_V2.R         # Estimación 3 CP — todos los modelos (versión final)
│   ├── runJAGS_exploratorio.R   # Script exploratorio inicial (modelo individual)
│   ├── runJAGS_ciclo_0cp.R      # Ciclo para 0 CP usando coda.samples
│   ├── runJAGS_1cp_exploratorio.R  # Versión exploratoria 1 CP
│   └── runJAGS_2cp_exploratorio.R  # Versión exploratoria 2 CP
│
├── 04_visualizacion/
│   ├── plot_lambdas.R              # Función de riesgo λ(t) — individual y comparativa
│   ├── plot_funcion_acumulada.R    # Función acumulada M(t) estimada vs observada
│   ├── plot_series_puntos_cambio.R # Series de tiempo con puntos de cambio marcados
│   ├── plot_sin_puntos_cambio.R    # Gráficas para el modelo sin puntos de cambio
│   └── plot_lambda_simple.R        # Versión simplificada de λ(t) (script inicial)
│
├── 05_diagnosticos_convergencia/
│   └── tests_convergencia.R     # Gelman-Rubin, Geweke, Heidel-Welch, Raftery-Lewis
│                                #   + gráficas de densidad, trazas, autocorrelación
│
├── 06_tablas_resultados/
│   └── tablas_resultados.R      # Tablas LaTeX: parámetros, DIC, WAIC, SDM,
│                                #   test de Ljung-Box por lag
│
└── 07_utilidades/
    └── crear_mcmc_list.R        # Convierte simulaciones de R2jags a mcmc.list
```

---

## Datasets analizados (10 modelos)

| # | Variable   | Ciudad   | Umbral (µg/m³ o ppb) |
|---|------------|----------|----------------------|
| 1 | Ozono      | Bogotá   | 51                   |
| 2 | Ozono      | Medellín | 51                   |
| 3 | PM10       | Bogotá   | 75 (Colombia)        |
| 4 | PM10       | Medellín | 75 (Colombia)        |
| 5 | PM10       | Bogotá   | 50 (OMS)             |
| 6 | PM10       | Medellín | 50 (OMS)             |
| 7 | PM2.5      | Bogotá   | 37 (Colombia)        |
| 8 | PM2.5      | Medellín | 37 (Colombia)        |
| 9 | PM2.5      | Bogotá   | 25 (OMS)             |
|10 | PM2.5      | Medellín | 25 (OMS)             |

---

## Flujo de trabajo recomendado

1. **Preparar datos** → `01_datos/`
2. **Estimar modelos** → `03_estimacion/` (usar versiones `_V2`)
   - Los resultados se guardan como `0cp.RData`, `1cp.RData`, `2cp.RData`, `3cp.RData`
3. **Visualizar resultados** → `04_visualizacion/`
4. **Revisar convergencia** → `05_diagnosticos_convergencia/`
5. **Generar tablas** → `06_tablas_resultados/`

> **Nota:** Cargar `07_utilidades/crear_mcmc_list.R` antes de los scripts de
> estimación y diagnóstico, ya que varios archivos dependen de esta función.

---

## Dependencias principales de R

```r
library(R2jags)      # Interfaz con JAGS
library(rjags)
library(coda)        # Objetos mcmc.list
library(bayesplot)   # Histogramas y gráficas MCMC
library(loo)         # WAIC y LOO
library(lattice)     # densityplot / xyplot para cadenas
library(xtable)      # Tablas LaTeX
library(readxl)      # Lectura datos Bogotá
library(dplyr)
library(writexl)
```
