## Procesamiento de los datos de Bogotá 
## Gráficas de las series de tiempo
## Cálculo de Excedencias 
## Aislamiento Preventivo Obligatorio

library(readxl)
data <- read_excel("rmcab_1h_ene 2018-dic 2020.xls", 
                   range = "A3:DB26309", col_names = FALSE)
data <- data[-3,] #Borrar fila de unidades
data[data == 'Sin Data'] <- NA #Pasar 'Sin Data' a NA

# PM2.5 ------------------------------------------------------------------------
pm25C <- data[ , grepl('PM2.5', data[2,])] #Seleccionar los datos de PM2.5 para 
                                           #todas las estaciones y todas las horas
pm25C <- cbind(data$...1, pm25C) # Añadir columna de fechas
colnames(pm25C) <- pm25C[1,]     # Cambiar nombres de las columnas
pm25C <- pm25C[-c(1,2),]         # Eliminar las filas 1 y 2 por no ser relevantes

# Formatear las fechas
pm25C$`Fecha & Hora` <- as.Date(pm25C$`Fecha & Hora`, format="%d/%m/%Y ")
colnames(pm25C)[1] <- 'Fecha'

# Pasar datos a números
pm25C[ , 2:19] <- apply(pm25C[ , 2:19], 2, function(x) as.numeric(as.character(x)))

# Promedios por día para cada estación
promedios <- list()
promedios <- aggregate(.~Fecha, pm25C, FUN=mean, na.rm = TRUE, na.action = na.pass)

# Todas las fechas posibles
fechaInicial <- as.Date("2018-01-01")
fechaFinal <- as.Date("2020-12-31")
fecha <- seq(from=fechaInicial, to=fechaFinal, by=1)

# Promedio de todas las estaciones para cada día
PM25 <- apply(promedios[,-1], MARGIN = 1, function(x) max(x,na.rm=T))
PM25 <- data.frame(fecha, PM25)


# Gráficas de las series de tiempo ---------------------------------------------

plot(as.Date(PM25$fecha), PM25$PM25,'l', xlab="Dates", ylab="PM2.5 [μg/m3]")
abline(h=25, col="blue")
abline(h=37, col="red")
abline(v=as.Date("2020-03-25"), col="green")
legend(x = "topleft", legend = c("WHO", "Colombia"), fill = c("blue", "red"), cex = 0.6)

# Excedencias ------------------------------------------------------------------

pm25who <- 25
pm25col <- 37

excePM25who_Bog <- which(PM25$PM25 > pm25who)
excePM25col_Bog <- which(PM25$PM25 > pm25col)

# Mandatory Preventive Isolation -----------------------------------------------

day <- sum(fecha < as.Date("2020-03-25"))
sum(excePM25who_Bog < day)
sum(excePM25who_Bog >= day)
sum(excePM25col_Bog < day)
sum(excePM25col_Bog >= day)
