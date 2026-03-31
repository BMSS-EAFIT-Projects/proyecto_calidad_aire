# Cumulative function

load("NHPP_3cp.RData")

plot_M_estimation <- function(n_modelo,
                               y_name, 
                              samples_n){
  
  d <- di[[n_modelo]]
  acum <- 0
  for(i in 1:1096){
    acum[i] <- sum(d<=i)
  }
  summary <- samples_n[[n_modelo]][["BUGSoutput"]][["summary"]]
  m <- summary[grep("m", row.names(summary)),1]
  m_2.5 <- summary[grep("m", row.names(summary)),3]
  m_97.5 <- summary[grep("m", row.names(summary)),7]
  
  plot(NA, type= 'n', xlab = "Days", ylab = y_name[n_modelo], lwd =2,  xlim = c(0,1096), 
       ylim = c(min(m_2.5),max(m_97.5)))
  
  lines(d,m, col="blue", lty = 'solid', lwd = 1)
  lines(d,m_2.5, col="red", lty = 'dashed')
  lines(d,m_97.5, col="red", lty = 'dashed')
  
  lines(acum, lty= 'solid', lwd = 1.5) 
  
  # Para d=4 PM 10 Col Medellín especial
  if(n_modelo == 4){
    plot(NA, type= 'n', xlab = "Days", ylab = y_name[n_modelo], lwd =2,  xlim = c(0,1096), 
         ylim = c(min(m_2.5),61))
    
    lines(c(d,1096),c(m,46), col="blue", lty = 'solid', lwd = 1)
    lines(c(d,1096),c(m_2.5, 34), col="red", lty = 'dashed')
    lines(c(d,1096),c(m_97.5,60), col="red", lty = 'dashed')
    
    lines(acum, lty= 'solid', lwd = 1.5)
  }
}



calculate_m <- function(n_modelo, samples_n = samples_resultados3) {
  #d <- di[[n_modelo]]
  d <- 1:1096
  su <- samples_n[[n_modelo]]$BUGSoutput$summary
  alpha <- su[grep("alpha", row.names(su)),][,1]
  beta <- su[grep("beta", row.names(su)),][,1]
  tau <- su[grep("tau", row.names(su)),][,1]
  
  t <- d
  # Validación de parámetros
  if(length(alpha) != length(beta)) {
    stop("alpha y beta deben tener la misma longitud")
  }
  
  puntos_cambio <- if(is.null(tau)) 0 else length(tau)
  n_segmentos <- puntos_cambio + 1
  
  if(length(alpha) != n_segmentos) {
    stop("El número de alphas no coincide con los puntos de cambio + 1")
  }
  if(length(beta) != n_segmentos) {
    stop("El número de betas no coincide con los puntos de cambio + 1")
  }
  
  # Función base m(t) para un segmento
  m_segment <- function(t, a, b) {
    (t/b)^a
  }
  
  # Caso sin puntos de cambio
  if(puntos_cambio == 0) {
    return(m_segment(t, alpha, beta))
  }
  
  # Inicializar vector de resultados
  result <- numeric(length(t))
  
  # Primer segmento (0 ≤ t < τ1)
  seg1 <- t < tau[1]
  result[seg1] <- m_segment(t[seg1], alpha[1], beta[1])
  
  # Segundo segmento (τ1 ≤ t < τ2 o t ≥ τ1 si solo hay 1 punto)
  if(puntos_cambio >= 1) {
    if(puntos_cambio == 1) {
      seg2 <- t >= tau[1]
    } else {
      seg2 <- t >= tau[1] & t < tau[2]
    }
    
    m_tau1_theta1 <- m_segment(tau[1], alpha[1], beta[1])
    result[seg2] <- m_tau1_theta1 + 
      m_segment(t[seg2], alpha[2], beta[2]) - 
      m_segment(tau[1], alpha[2], beta[2])
  }
  
  # Tercer segmento (solo si hay al menos 2 puntos)
  if(puntos_cambio >= 2) {
    if(puntos_cambio == 2) {
      seg3 <- t >= tau[2]
    } else {
      seg3 <- t >= tau[2] & t < tau[3]
    }
    
    m_tau1_theta1 <- m_segment(tau[1], alpha[1], beta[1])
    m_tau1_theta2 <- m_segment(tau[1], alpha[2], beta[2])
    m_tau2_theta2 <- m_segment(tau[2], alpha[2], beta[2])
    
    result[seg3] <- m_segment(t[seg3], alpha[3], beta[3]) - 
      m_segment(tau[2], alpha[3], beta[3]) +
      (m_tau2_theta2 - m_tau1_theta2) +
      m_tau1_theta1
  }
  
  # Cuarto segmento (solo si hay 3 puntos)
  if(puntos_cambio == 3) {
    seg4 <- t >= tau[3]
    
    m_tau1_theta1 <- m_segment(tau[1], alpha[1], beta[1])
    m_tau1_theta2 <- m_segment(tau[1], alpha[2], beta[2])
    m_tau2_theta2 <- m_segment(tau[2], alpha[2], beta[2])
    m_tau2_theta3 <- m_segment(tau[2], alpha[3], beta[3])
    m_tau3_theta3 <- m_segment(tau[3], alpha[3], beta[3])
    
    result[seg4] <- m_segment(t[seg4], alpha[4], beta[4]) - 
      m_segment(tau[3], alpha[4], beta[4]) +
      (m_tau3_theta3 - m_tau2_theta3) +
      (m_tau2_theta2 - m_tau1_theta2) +
      m_tau1_theta1
  }
  
  return(result)
}


# Pruebas----------------------

bog <- 1
med <- 4

y_name <- c("PM 2.5 Colombia 37", 
            "PM 2.5 Colombia 37", "PM2.5 WHO 25","PM2.5 WHO 25")
par(mfcol=c(1,1))
plot_M_estimation(bog, y_name, samples_resultados3)
lines(calculate_m(1), lwd = 2, col = "magenta")


# Guardar gráfica --------------------------
name_img <- "esti_Med_PM2.5WHO_1CP.png"
dir <- "/Users/karen/Desktop/analisisDatos/Graficas/Graficas_3_puntos/"
path <- paste(dir, name_img, sep="")
png(filename = path, width = 1700, height = 1200, res=280)
#par(cex = 1)
par(mar = c(5.1, 4.1, 2.1, 2.1))
#par(oma = c(0, 0, 0, 0))
plot_M_estimation(10, y_name, samples_resultados1)
dev.off()



