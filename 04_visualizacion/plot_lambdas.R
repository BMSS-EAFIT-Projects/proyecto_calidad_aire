# Para 3 puntos de cambio solamente y una sola ciudad
plot_lambda_single <- function(n_modelo,
                               y_name = "Lambda",
                               city_name = "Ciudad",
                               line_type = "solid",
                               point_color = "black",
                               line_width = 1.5) {
  
  #d = di[[n_modelo]]
  d=1:1096
    
  su <- samples_resultados3[[n_modelo]]$BUGSoutput$summary
  
  alpha <- su[grep("alpha", row.names(su)),][,1]
  beta <- su[grep("beta", row.names(su)),][,1]
  tau <- su[grep("tau", row.names(su)),][,1]
  
  
  # ---- Calcular lambda ----
  lambda <- ifelse(d < tau[1], 
                   (alpha[1]/beta[1]) * (d/beta[1])^(alpha[1]-1),
                   ifelse(d < tau[2], 
                          (alpha[2]/beta[2]) * (d/beta[2])^(alpha[2]-1),
                          ifelse(d < tau[3], 
                                 (alpha[3]/beta[3]) * (d/beta[3])^(alpha[3]-1),
                                 (alpha[4]/beta[4]) * (d/beta[4])^(alpha[4]-1))))
  
  # ---- Configurar gráfica ----
  plot(NA, xlim = range(d), ylim = range(lambda), 
       xlab = "Days", ylab = y_name, type = "n")
  
  # ---- Graficar segmentos ----
  sep <- sapply(tau, function(t) which(d >= t)[1])
  sep <- c(1, sep, length(d)+1)
  
  for(i in 1:(length(sep)-1)) {
    idx <- seq(sep[i], sep[i+1]-1)
    lines(d[idx], lambda[idx], col = point_color, lty = line_type, lwd = line_width)
  }
  
  # ---- Añadir puntos de cambio ----
  # abline(v = tau, col = adjustcolor(point_color, alpha.f = 0.4), lty = 3)
  
}

y_name <- c("O3_51", "O3_51", "PM10_75", "PM10_75", "PM10_50", "PM10_50",
            "PM2.5_37", "PM2.5_37", "PM2.5_25","PM2.5_25")
bog <- 3
med <- 4

# Guardar lambda Bog
name_img = "lambda_Bog_PM2.5WHO_3CP.png"
dir = "/Users/karen/Desktop/analisisDatos/Graficas/Graficas_3_puntos/"
path = paste(dir, name_img, sep="")

png(filename=path, width = 614, height = 419)
plot_lambda_single(bog, y_name = y_name[bog], line_type = "dashed", line_width = 2)
dev.off()

# Guardar lambda Med
name_img = "lambda_Med_PM2.5WHO_3CP.png"
dir = "/Users/karen/Desktop/analisisDatos/Graficas/Graficas_3_puntos/"
path = paste(dir, name_img, sep="")

png(filename=path, width = 614, height = 419)
plot_lambda_single(med, y_name = y_name[med], line_type = "solid", line_width = 1.5)
dev.off()

#Lambda Bogota
plot_lambda_single(bog, y_name = y_name, line_type = "dashed", line_width = 2)

#Lambda Med
plot_lambda_single(med, y_name = y_name, line_type = "solid", line_width = 1.5)





# Comparación entre las dos ciudades ---------------

# Esta función calcula una gráfica comparativa de los lambdas para ambas
# ciudades, requiere el modelo (dataset del 1 al 10), el num de puntos
# de cambio (0 a 3) para ciudad
plot_lambda_comparison <- function(modelo_bog, modelo_med, 
                                   puntos_cambio_bog, 
                                   puntos_cambio_med,
                                   y_label = "Lambda",
                                   ciudad_bog = "Bogotá",
                                   ciudad_med = "Medellín",
                                   colors = c("red", "blue")) {
  
  samples_list <- list(samples_resultados0, samples_resultados1, samples_resultados2, samples_resultados3)
  
  # Función para extraer tau
  extraer_tau <- function(su, puntos_cambio) {
    if(puntos_cambio == 0) return(NULL)
    
    tau_rows <- grep("tau", row.names(su))
    if(length(tau_rows) == 0) return(NULL)
    
    # Manejo especial para 1 punto de cambio
    if(puntos_cambio == 1) {
      tau <- su[tau_rows, "mean"]
      if(is.matrix(tau)) tau <- tau[,1]  # Si es matriz, extraer primera columna
      return(c(tau))  # Asegurar que es vector
    } else {
      return(su[tau_rows, "mean"])
    }
  }
  
  # Extraer parámetros para Bogotá
  su_bog <- samples_list[[puntos_cambio_bog+1]][[modelo_bog]]$BUGSoutput$summary
  
  alpha_bog <- su_bog[grep("alpha", row.names(su_bog)), "mean"]
  beta_bog <- su_bog[grep("beta", row.names(su_bog)), "mean"]
  tau_bog <- extraer_tau(su_bog, puntos_cambio_bog)
  
  # Extraer parámetros para Medellín
  su_med <- samples_list[[puntos_cambio_med+1]][[modelo_med]]$BUGSoutput$summary
  
  alpha_med <- su_med[grep("alpha", row.names(su_med)), "mean"]
  beta_med <- su_med[grep("beta", row.names(su_med)), "mean"]
  tau_med <- extraer_tau(su_med, puntos_cambio_med)
  
  # Días
  d <- 1:1096
  
  # Función para calcular lambda mejorada
  calcular_lambda <- function(d, alpha, beta, tau, puntos_cambio) {
    if(puntos_cambio == 0) {
      return((alpha/beta) * (d/beta)^(alpha-1))
    } else {
      lambda <- numeric(length(d))
      
      # Primer segmento (0 <= t < tau1)
      seg1 <- d < tau[1]
      lambda[seg1] <- (alpha[1]/beta[1]) * (d[seg1]/beta[1])^(alpha[1]-1)
      
      # Segmentos intermedios (solo si hay más de 1 punto)
      if(puntos_cambio > 1) {
        for(j in 2:puntos_cambio) {
          seg <- d >= tau[j-1] & d < tau[j]
          lambda[seg] <- (alpha[j]/beta[j]) * (d[seg]/beta[j])^(alpha[j]-1)
        }
      }
      
      # Último segmento (tau_n <= t)
      seg_last <- d >= tau[puntos_cambio]
      lambda[seg_last] <- (alpha[puntos_cambio+1]/beta[puntos_cambio+1]) * 
        (d[seg_last]/beta[puntos_cambio+1])^(alpha[puntos_cambio+1]-1)
      
      return(lambda)
    }
  }
  
  # Calcular lambdas
  lambda_bog <- calcular_lambda(d, alpha_bog, beta_bog, tau_bog, puntos_cambio_bog)
  lambda_med <- calcular_lambda(d, alpha_med, beta_med, tau_med, puntos_cambio_med)
  
  # Configurar gráfica
  ylim <- range(c(lambda_bog, lambda_med), na.rm = TRUE)
  
  plot(NA, xlim = range(d), ylim = ylim, 
       xlab = "Days", ylab = y_label, type = "n")
  
  # Función para graficar segmentos
  graficar_segmentos <- function(d, lambda, tau, col, lty, lwd) {
    if(!is.null(tau) && length(tau) > 0) {
      sep <- sapply(tau, function(t) which(d >= t)[1])
      sep <- c(1, sep, length(d)+1)
      
      for(i in 1:(length(sep)-1)) {
        idx <- seq(sep[i], sep[i+1]-1)
        lines(d[idx], lambda[idx], col = col, lty = lty, lwd = lwd)
      }
    } else {
      lines(d, lambda, col = col, lty = lty, lwd = lwd)
    }
  }
  
  # Graficar Bogotá
  graficar_segmentos(d, lambda_bog, tau_bog, 
                     col = colors[1], lty = "dashed", lwd = 2)
  
  # Graficar Medellín
  graficar_segmentos(d, lambda_med, tau_med, 
                     col = colors[2], lty = "solid", lwd = 1.5)
  
  # Añadir puntos de cambio
  if(!is.null(tau_bog)) abline(v = tau_bog, col = adjustcolor(colors[1], alpha.f = 0.3), lty = 3)
  if(!is.null(tau_med)) abline(v = tau_med, col = adjustcolor(colors[2], alpha.f = 0.3), lty = 3)
  
  # Leyenda
  legend("topright", legend = c(ciudad_bog, ciudad_med),
         col = colors, lty = c("dashed", "solid"), 
         lwd = c(2, 1.5), cex = 0.6)
}

y_name <- c("Ozone 51", "Ozone 51", "PM10 Colombia 75", "PM10 Colombia 75", 
            "PM10 WHO 50", "PM10 WHO 50","PM 2.5 Colombia 37", 
            "PM 2.5 Colombia 37", "PM2.5 WHO 25","PM2.5 WHO 25")



# Graficar comparación 
plot_lambda_comparison(
  modelo_bog = bog,
  modelo_med = med,
  puntos_cambio_bog = puntos_bog,
  puntos_cambio_med = puntos_med,
  y_label = y_name[bog],
  ciudad_bog = "Bogotá",
  ciudad_med = "Medellín"
)

# Índices de modelos 
bog <- 9   # Modelo para Bogotá
puntos_bog <- 1
med <- 10 # Modelo para Medellín
puntos_med <- 1

# Guardar gráfica 
name_img <- "final_lambda_PM25who.png"
dir <- "/Users/karen/Desktop/analisisDatos/Graficas/Finales/"
path <- paste(dir, name_img, sep="")

png(filename = path, width = 1700, height = 1200, res = 280)
par(mar = c(5.1, 4.1, 2.1, 2.1))
plot_lambda_comparison(
  modelo_bog = bog,
  modelo_med = med,
  puntos_cambio_bog = puntos_bog,
  puntos_cambio_med = puntos_med,
  y_label = y_name[bog],
  ciudad_bog = "Bogotá",
  ciudad_med = "Medellín"
)
dev.off()



