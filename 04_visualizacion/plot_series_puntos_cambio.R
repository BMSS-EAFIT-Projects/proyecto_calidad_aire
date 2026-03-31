
obtener_tau <- function(modelo, puntos){
  samples_list <- list(samples_resultados0, samples_resultados1, samples_resultados2, samples_resultados3)
  
  su <- samples_list[[puntos+1]][[modelo]]$BUGSoutput$summary
  if(puntos == 1){
    tau <- su[grep("tau", row.names(su)),][1]
  } else{
    tau <- su[grep("tau", row.names(su)),][,1]
  }
  tau <- round(tau,0)
  return(tau)
}

modelo <- 9
puntos_cambio <- 1

taus <- obtener_tau(modelo,puntos_cambio)
taus
fecha[taus]

# BOGOTA
plot(as.Date(OZONO$fecha), OZONO$OZONO,'l', xlab="Dates", ylab="Ozone [ppb]")
abline(h=51, col="blue", lwd = 1.5)
abline(h=51, col="red", lwd = 1.5) 
abline(v=as.Date("2020-03-25"), col="green", lwd = 2)
abline(v=as.Date(fecha[taus]), col = "darkorchid", lwd = 2)
legend(x = "topleft", legend = c("WHO", "Colombia", "Change Points"), fill = c("blue", "red", "darkorchid"), cex = 0.6)

plot(as.Date(PM10$fecha), PM10$PM10,'l', xlab="Dates", ylab="PM10 [μg/m3]")
abline(h=50, col="blue", lwd = 1.5)
abline(h=75, col="red", lwd = 1.5)
abline(v=as.Date("2020-03-25"), col="green")
abline(v=as.Date(fecha[taus]), col = "orange", lwd = 2)
legend(x = "topleft", legend = c("WHO", "Colombia"), fill = c("blue", "red"), cex = 0.6)

plot(as.Date(PM25$fecha), PM25$PM25,'l', xlab="Dates", ylab="PM2.5 [μg/m3]")
abline(h=25, col="blue", lwd = 1.5)
abline(h=37, col="red", lwd = 1.5)
abline(v=as.Date("2020-03-25"), col="green")
abline(v=as.Date(fecha[taus]), col = "orange", lwd = 2)
legend(x = "topleft", legend = c("WHO", "Colombia"), fill = c("blue", "red"), cex = 0.6)

# MEDELLIN
plot(as.Date(ozone_Med$fecha), ozone_Med$Ozono,'l', xlab="Dates", ylab="Ozone [ppb]")
abline(h=51, col="blue", lwd = 1.5)
abline(h=51, col="red", lwd = 1.5)
abline(v=as.Date("2020-03-25"), col="green")
abline(v=as.Date(fecha[taus]), col = "orange", lwd = 2)
legend(x = "topleft", legend = c("WHO", "Colombia"), fill = c("blue", "red"), cex = 0.6)

plot(as.Date(PM10_Med$fecha), PM10_Med$PM10,'l', xlab="Dates", ylab="PM10 [μg/m3]")
abline(h=50, col="blue", lwd = 1.5)
abline(h=75, col="red", lwd = 1.5)
abline(v=as.Date("2020-03-25"), col="green")
abline(v=as.Date(fecha[taus]), col = "orange", lwd = 2)
legend(x = "topleft", legend = c("WHO", "Colombia"), fill = c("blue", "red"), cex = 0.6)

plot(as.Date(PM25_Med$fecha), PM25_Med$PM25,'l', xlab="Dates", ylab="PM2.5 [μg/m3]")
abline(h=25, col="blue", lwd = 1.5)
abline(h=37, col="red", lwd = 1.5)
abline(v=as.Date("2020-03-25"), col="green")
abline(v=as.Date(fecha[taus]), col = "orange", lwd = 2)
legend(x = "topleft", legend = c("WHO", "Colombia"), fill = c("blue", "red"), cex = 0.6)

# Guardar gráfica 
d_names <- c("Ozone_Bog", "Ozone_Med", "PM10col_Bog", "PM10col_Med", "PM10who_Bog",
             "PM10who_Med", "PM25col_Bog", "PM25col_Med", "PM25who_Bog", "PM25who_Med")


# Correr los puntos de cambio y el modelo-----------
mod_col <- 8
pc_col <- 3

mod_who <- 10
pc_who <- 1


name_img <- paste(d_names[mod_col], "_", d_names[mod_who], ".png", sep="")
dir <- "/Users/karen/Desktop/analisisDatos/Graficas/Series_con_puntos/"
path <- paste(dir, name_img, sep="")
png(filename = path, width = 1800, height = 1200, res = 280)
#c(bottom, left, top, right)
par(mar = c(5.1, 4.1, 2.1, 2.1))

# Cambiar gráfica
plot(as.Date(PM25_Med$fecha), PM25_Med$PM25,'l', xlab="Dates", ylab="PM2.5 [μg/m3]")
abline(h=25, col="blue", lwd = 1.5)
abline(h=37, col="red", lwd = 1.5)

abline(v=as.Date("2020-03-25"), col="green2", lwd = 1.5)
taus <- obtener_tau(mod_col,pc_col)
abline(v=as.Date(fecha[taus]), col = "darkorange", lwd = 2, lty = 2)

taus <- obtener_tau(mod_who,pc_who)
abline(v=as.Date(fecha[taus]), col = "darkorchid", lwd = 2, lty = 2)

dev.off()


# Modified legend command
legend(x = "right", 
       legend = c("WHO Threshold", "Colombia Threshold", "Mandatory Preventive Isolation", "Change Points WHO", "Change Points Col"), 
       col = c("blue", "red",  "green2","darkorchid", "darkorange" ), 
       lwd = 2, 
       cex = 0.6,
       xpd = TRUE,  # Allows drawing outside plot area
       inset = c(-0.5, 0),
       lty = c(1,1,1,3,3)) 

name_img <- "leyenda_horizontal.png"
dir <- "/Users/karen/Desktop/analisisDatos/Graficas/Series_con_puntos/"
path <- paste(dir, name_img, sep="")
# Crear un dispositivo gráfico PNG para la leyenda
png(path, width = 3200, height = 180, res=200)

# Configurar un gráfico vacío
par(mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

# Crear la leyenda horizontal
legend("center", 
       legend = c("WHO Threshold", "Colombia Threshold", 
                  "Change Points WHO", "Change Points Col", "Mandatory Preventive Isolation"), 
       col = c("blue", "red", "darkorchid", "darkorange", "green2"), 
       lwd = 2, 
       cex = 1.1,
       bty = "n",  # Sin caja alrededor
       horiz = TRUE,  # Leyenda horizontal
       lty = c(1, 1, 1, 3, 3)
       )

# Cerrar el dispositivo gráfico
dev.off()
