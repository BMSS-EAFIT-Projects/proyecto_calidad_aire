
#### Cumulative Function
bog <- 9
med <- 10
y_name <- "PM2.5_25"

d_M <- di[[med]]
d_B <- di[[bog]]

acum_M <- 0
for(i in 1:1096){
  acum_M[i] <- sum(d_M<=i)
}

acum_B <- 0
for(i in 1:1096){
  acum_B[i] <- sum(d_B<=i)
}

m_M <- summary_models[[med]]$statistics[1:length(d_M)+2,1]
m_B <- summary_models[[bog]]$statistics[1:length(d_B)+2,1]

m_M1 <- summary_models[[med]]$quantiles[1:length(d_M)+2,1]
m_M2 <- summary_models[[med]]$quantiles[1:length(d_M)+2,5]

m_B1 <- summary_models[[bog]]$quantiles[1:length(d_B)+2,1]
m_B2 <- summary_models[[bog]]$quantiles[1:length(d_B)+2,5]

plot(0, type= 'l', xlab = "Days", ylab = y_name, lwd =2,  xlim = c(0,1096), 
     ylim = c(min(m_M1),max(m_M2)))

plot(0, type= 'l', xlab = "Days", ylab = y_name, lwd =2,  xlim = c(0,1096), 
     ylim = c(min(c(acum_B,acum_M)), max(c(acum_B,acum_M))))

lines(d_M, acum_M, lty= 'solid', lwd =1.5) 
lines(acum_M, lty = 'dashed', lwd = 2)

legend(x = "topleft", legend = c("Medellin", "Bogota"), lty=c("solid","dashed"), 
       cex = 0.6)

lines(d_M,m_M, col="blue", lty = 'solid', lwd = 1)
lines(d_M,m_M1, col="red", lty = 'dashed')
lines(d_M,m_M2, col="red", lty = 'dashed')

lines(d_B,m_B, col="blue", lty = 'solid', lwd = 1)
lines(d_B,m_B1, col="red", lty = 'dashed')
lines(d_B,m_B2, col="red", lty = 'dashed')


### Special plot for di[[4]]
plot(0, type= 'l', xlab = "Days", ylab = y_name, lwd =2,  xlim = c(1,1096), 
     ylim = c(0,60))
lines(c(1,d_M,1096), c(0.7,m_M, 44.83102), col="blue", lty = 'solid', lwd = 1)
lines(c(1,d_M,1096), c(-0.5,m_M1, 33), col="red", lty = 'dashed')
lines(c(1,d_M,1096), c(2.5,m_M2, 58.65082), col="red", lty = 'dashed')


# Parameters para LAMBDA

alpha_M <- summary_models[[med]]$statistics[1,1]
beta_M <- summary_models[[med]]$statistics[2,1]

alpha_B <- summary_models[[bog]]$statistics[1,1]
beta_B <- summary_models[[bog]]$statistics[2,1]


#### Lambda Rate function 

lambda_M <- 0
for(i in 1:length(d_M)){
  lambda_M[i] <- (alpha_M/beta_M)*(d_M[i]/beta_M)^(alpha_M-1)
}


lambda_B <- 0
for(i in 1:length(d_B)){
  lambda_B[i] <- (alpha_B/beta_B)*(d_B[i]/beta_B)^(alpha_B-1)
}

# Para grafica individual (cambiar limites de y)
plot(0, type= 'l', xlab = "Days", ylab = y_name, lwd =2,  xlim = c(0,1096), 
     ylim = c(min(lambda_M),max(lambda_M)))

# Para gráfica comparativa
plot(0, type= 'l', xlab = "Days", ylab = y_name, lwd =2,  xlim = c(0,1096), 
     ylim = c(min(c(lambda_B,lambda_M)), max(c(lambda_B,lambda_M))))

lines(d_M, lambda_M, lty = 'solid', lwd = 1.5)
lines(d_B, lambda_B, lty = 'dashed', lwd =2)

legend(x = "topright", legend = c("Medellin", "Bogota"), lty=c("solid","dashed"), 
       cex = 0.6)

# Special plot for di[[4]]
special_d_M <- c(1,2,3,4,5,6,7,8,9,10,11,12,15,20,25,30,35,40,45, d_M, 1000,1096)
special_lambda_M <- 0
for(i in 1:length(special_d_M)){
  special_lambda_M[i] <- (alpha_M/beta_M)*(special_d_M[i]/beta_M)^(alpha_M-1)
}


lines(special_d_M,special_lambda_M, lty = 'solid', lwd = 1.5)



