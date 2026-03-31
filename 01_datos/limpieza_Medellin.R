## Medellin Data

# Read Medellin data
PM25_Med <- read.csv("datosMedellin/PM25_resumen.csv")

PM25_Med <- subset(PM25_Med, (fecha >= "2018-01-01") & (fecha <= "2020-12-31"))

rownames(PM25_Med) <- NULL

# Exceedances ------------------------------------------------------------------

pm25who <- 25
pm25col <- 37

excePM25who_Med <- which(PM25_Med$PM25 > pm25who)
excePM25col_Med <- which(PM25_Med$PM25 > pm25col)

# Mandatory Preventive Isolation -----------------------------------------------

fecha <- as.Date(PM25_Med$fecha)
day <- sum(fecha < as.Date("2020-03-25"))
sum(excePM25who_Med < day)
sum(excePM25who_Med >= day)
sum(excePM25col_Med < day)
sum(excePM25col_Med >= day)

# Gráficas de las series de tiempo ---------------------------------------------

plot(as.Date(PM25_Med$fecha), PM25_Med$PM25,'l', xlab="Dates", ylab="PM2.5 [μg/m3]")
abline(h=25, col="blue")
abline(h=37, col="red")
abline(v=as.Date("2020-03-25"), col="green")
legend(x = "topleft", legend = c("WHO", "Colombia"), fill = c("blue", "red"), cex = 0.6)
