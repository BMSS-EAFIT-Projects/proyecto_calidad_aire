#Convertir a objeto mcmc.list
library(coda)  # Asegúrate de cargar la librería coda

crear_mcmc_list <- function(output_cadenas, n_chains) {
  # Verificar que el número de cadenas sea válido
  if (n_chains <= 0) {
    stop("El número de cadenas (n_chains) debe ser mayor que 0.")
  }
  
  # Calcular el tamaño de cada cadena
  total_iteraciones <- nrow(output_cadenas)
  iteraciones_por_cadena <- total_iteraciones / n_chains
  
  # Verificar que el número de iteraciones sea divisible entre el número de cadenas
  if (iteraciones_por_cadena != round(iteraciones_por_cadena)) {
    stop("El número de iteraciones no es divisible entre el número de cadenas.")
  }
  
  # Crear una lista para almacenar las cadenas
  lista_cadenas <- list()
  
  # Dividir output_cadenas en n_chains partes y convertirlas en objetos mcmc
  for (i in 1:n_chains) {
    inicio <- (i - 1) * iteraciones_por_cadena + 1
    fin <- i * iteraciones_por_cadena
    cadena <- mcmc(window(output_cadenas, inicio, fin))
    lista_cadenas[[i]] <- cadena
  }
  
  # Convertir la lista en un objeto mcmc.list
  mcmc_list <- mcmc.list(lista_cadenas)
  
  return(mcmc_list)
}
