# Lambda Rate function 
# d1 Med d2 Bog

d1 <- excePM10col_Med
d2 <- excePM25who_Bog

alpha1 <- summary_models[[4]]$statistics[1,1]
beta1 <- summary_models[[4]]$statistics[2,1]

lambd1 <- 0
for(i in 1:length(d1)){
  lambd1[i] <- (alpha1/beta1)*(d1[i]/beta1)^(alpha1-1)
}

alpha2 <- summary_models[[9]]$statistics[1,1]
beta2 <- summary_models[[9]]$statistics[2,1]

lambd2 <- 0
for(i in 1:length(d2)){
  lambd2[i] <- (alpha2/beta2)*(d2[i]/beta2)^(alpha2-1)
}

y_name <- "PM2.5_25"

plot(d1,lambd1, type= 'p', xlab = "Days", ylab = y_name, xlim = c(0,1096), ylim=c(0,0.2)) 
plot(d2,lambd2, type= 'l', lty = "dashed", xlab = "Days", ylab = y_name, lwd = 1.5) 


plot(d1,lambd1, type= 'l', xlab = "Days", ylab = y_name, 
     ylim=c(min(c(lambd1,lambd2)),max(c(lambd1,lambd2))), xlim=c(0,1096)) 
lines(d2,lambd2, lty = 'dashed', lwd = 1)
legend(100,1, legend = c("Medellin", "Bogota"), lty=c("solid","dashed"), cex = 0.6)
