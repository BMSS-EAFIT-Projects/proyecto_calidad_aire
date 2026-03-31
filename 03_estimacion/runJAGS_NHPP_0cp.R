# Estimation zero change point

library('xtable')
library(rjags)
library(lattice)
library(stats)
library(bayesplot)
library(loo)
library(R2jags)

options(mc.cores = 3)

set.seed(42) 


t <- 1096 # Total days
di <- list(excePM25col_Bog, excePM25col_Med, 
           excePM25who_Bog, excePM25who_Med)
d_names <- c("PM25col_Bog", "PM25col_Med", "PM25who_Bog", "PM25who_Med")

samples_resultados0 <- list()


t1 <- Sys.time()
for(i in 1:4){
  d <-  di[[i]]
  niter <- 100000
  nthin <- 10
  
  line_data <- list("K" = length(d), "T"= t, "d" = d)
  samples_resultados0[[i]] <- jags(data = line_data,
                                   parameters.to.save = c("alpha", "beta", "m", "phi2"),
                                   n.chains = 2,
                                   n.iter = niter,
                                   n.burnin = 20000,
                                   n.thin = nthin,
                                   model.file = "02_modelos_JAGS/bug_files/NHPP_0cp.bug")
  
  print(paste("Listo modelo ", i))
  save(samples_resultados0, file="NHPP_0cp.RData")
}

Sys.time() - t1

load("NHPP_0cp.RData")


mod = 1
# Parameters
for(j in mod){
  d <-  di[[j]]
  samples <- samples_resultados0[[j]]$BUGSoutput
  y_name <- d_names[mod]
}


# Gráfica para modelo actual --------------------------
# Acumulada de los datos
acum <- 0
for(i in 1:1096){
  acum[i] <- sum(d<=i)
}  
# Acumulada estimada
su <- samples$summary
# mean de la m
mm <- su[grep("m", row.names(su)),1]

plot(0, type= 'l', xlab = "Days", ylab = y_name, lwd =2,  xlim = c(0,1096), ylim = c(0,max(c(mm,acum))))
lines(d,mm)
lines(acum, lty = 'dashed', lwd = 2)

sum(abs(acum[d]-mm))

# Get mean, SD and confidence intervals of alpha, beta and tau

su[grep("alpha", row.names(su)),]
su[grep("beta", row.names(su)),]


# Convergence Tests ------------------------------------------------------------

#Obtener las simulaciones
alphas <- samples$sims.list$alpha
colnames(alphas)<-c("alpha1")
betas <- samples$sims.list$beta
colnames(betas)<-c("beta1")


# Convertir a objeto mcmc.list
alphas <- crear_mcmc_list(alphas, samples$n.chains)
betas <- crear_mcmc_list(betas, samples$n.chains)


## Histograma
mcmc_hist(alphas)
mcmc_hist(betas)


## Alpha
trellis.par.set(strip.background=list(col="white"))

densityplot(alphas, ylab = '', main = 'Density of alpha')
xyplot(alphas, main = 'Trace of alpha')
autocorr.plot(alphas)
r
gelman.plot(alphas)
gelman.plot(alphas, autoburnin = F)
geweke.plot(alphas)
r
raftery.diag(alphas)
heidel.diag(alphas)
gelman.diag(alphas)[1]
geweke.diag(alphas)

## Beta
densityplot(betas, ylab = '', main = 'Density of beta')
xyplot(betas, main = 'Trace of beta')
autocorr.plot(betas)
r
gelman.plot(betas)
gelman.plot(betas, autoburnin = F)
geweke.plot(betas)
r
raftery.diag(betas)
heidel.diag(betas)
gelman.diag(betas)[1]
geweke.diag(betas)



# Obtener DIC ------------------------------------------------------------------
samples$DIC

# Obtener WAIC
loglik <-samples$sims.list$phi2
waic0 <- waic(loglik)
waic0
loo(loglik)



