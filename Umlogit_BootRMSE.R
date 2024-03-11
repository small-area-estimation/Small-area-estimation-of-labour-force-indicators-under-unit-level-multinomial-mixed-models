



## This R code corresponds to the section "Application to real data" of the paper 
## Small area estimation of labour force indicators under unit-level multinomial mixed models.
## This is the script to obtain the bootstrap estimation of rmse for the Multinomial Mixed Model.


## AUTHORS: Bugallo M., Esteban M.D., Hozba T., Morales D., Perez A.



rm(list=ls())
set.seed(2509)

source("function_MLE_Laplace.R")
source('AuxFun.R')

library(car)
library(nnet)
library(mefa)
library(lme4)
library(mvtnorm)
library(cubature)
library(nloptr)
library(fastDummies)


## A unit-level multinomial logit mixed model by province:sex
# Bootstrap estimation of the MSE

## EPA_2021T1
EPA_2021T1 <- read.csv("EPA_2021T1.csv", sep=';', header=T) 
EPA_2021T1_new <- EPA_2021T1[which(EPA_2021T1$EDAD1>15), ] 
attach(EPA_2021T1_new)

# Target vector variable
lab <- recode(AOI, " c('3','4') = 'employed'; c('5','6') = 'unemployed'; c('7','8','9') = 'inactive' ", 
              as.factor=TRUE, levels = c('employed','unemployed', 'inactive'))
lab1 <- as.numeric(lab=='employed')
lab2 <- as.numeric(lab=='unemployed')
lab3 <- as.numeric(lab=='inactive')

# Auxiliary variables
age4 <- cut(EDAD1, breaks=c(0,45,55,64,Inf), include.lowest=TRUE, labels=c("1","2","3","4"))
X <- model.matrix(multinom(lab ~ age4))

# Census
Censo_2021T1 <- read.csv('Censo_2021T1.csv', sep=',', header=T)

# Domains and sample sizes
n <- length(lab1)
D <- length(unique(PROV))*length(unique(SEXO1))

nd <- data.frame(aggregate(rep(1, n), by=list(PROV, SEXO1), sum))$x
ndk <- data.frame(aggregate(rep(1, n), by=list(PROV, SEXO1, age4), sum))$x

Ndk <- Censo_2021T1$Conteo + ndk

# Fitting ML algorithm
opts <- list( "algorithm" = "NLOPT_LN_NEWUOA_BOUND", maxeval=1000, print_level=0)

# Generation of the domain index d
dd <- rep(NA, n)
for (i in 1:n){
  if (SEXO1[i]==1){
    dd[i] <- PROV[i]
  } else
    dd[i] <- PROV[i] + 52
}

# Explanatory variables: age4 (lab1 and lab2)
matriz <- data.frame(X, dd)
rownames(matriz) <- NULL

X <- list()
for(d in 1:D) {
  X[[d]] <- vector()
  
  for(j in 1:nd[d]) {
    X[[d]] <- rbind(X[[d]], t(c(1, age4[matriz$dd==d][j])))
  }
  colnames(X[[d]]) <- c('Intercept', 'age4')
  # One hot encoder
  X[[d]] <- dummy_cols(X[[d]], select_columns = c("age4"))[,c(1,4,5,6)] # age4_1 as reference
  X[[d]] <- sapply(X[[d]], 'as.numeric')
  colnames(X[[d]]) <- NULL # Intecerpt, age4_2, age4_3, age4_4
}

x1 <- x2 <- X

# Target variables
ydj <- ydj2 <- list()
for(d in 1:D) {
  ydj[[d]] <- vector()
  
  for(j in 1:nd[d]) {
    ydj[[d]] <- rbind(ydj[[d]], t(c(lab1[matriz$dd==d][j], lab2[matriz$dd==d][j], lab3[matriz$dd==d][j])))
  }
  ydj2[[d]] <- ydj[[d]][,c(1,2)]
}

nudj <- 1
beta1.0 <- c(0.81, 0.36, -0.98, 4.44) # Obtained fitting a multinomial model without random effects (see Umlogit.R script)
beta2.0 <- c(-0.78, -0.01, -1.40, -5.94) # Obtained fitting a multinomial model without random effects (see Umlogit.R script)
beta.0 <- c(beta1.0, beta2.0)

p1 <- length(beta1.0)
p2 <- length(beta2.0)
p <- p1 + p2

phi1.0 <- 1
phi2.0 <- 1
phi.vect <- c(phi1.0, phi2.0)
phi <- as.matrix(bdiag(phi1.0, phi2.0))
q <- nrow(phi) + 1

seeds <- param.0 <- c(beta.0, phi1.0, phi2.0)

nud <- list()
for(d in 1:D) {
  nud[[d]] <- rep(1, nd[d])
}


# ML-Laplace aproximation

seeds.u <- split(matrix(0, D, 2), 1:D)  
opts <- list( "algorithm" = "NLOPT_LN_NEWUOA_BOUND", maxeval=1000, print_level=0)   
mle.L <- MLE_Laplace(x1, x2, ydj, nud, seeds, seeds.u, max.iter=1000, epsilon.u=0.001, epsilon.par=0.0001)   

param.est <- as.matrix(mle.L$par.est)
udk.pred <- mle.L$udk.pred

udk1.pred <- udk2.pred<-rep(NA, D)
for(d in 1:D) {
  udk1.pred[d] <- udk.pred[[d]][1]
  udk2.pred[d] <- udk.pred[[d]][2]
}

### Plug-in predictions of averages by domains
ddf <- as.matrix(data.frame(rep(1:52, 8), 
                            rep(c(rep(1,52),rep(6,52)),4), rep(1:104,4),
                            rbind(do.call("rbind", replicate(D, c(1,0,0,0), simplify = F)),
                                  do.call("rbind", replicate(D, c(1,1,0,0), simplify = F)),
                                  do.call("rbind", replicate(D, c(1,0,1,0), simplify = F)),
                                  do.call("rbind", replicate(D, c(1,0,0,1), simplify = F)))))
colnames(ddf) <- c('prov', 'sex', 'dominio', 'intercept', 'age42', 'age43', 'age44')

pdj_employed <- exp(ddf[, 4:7]%*%param.est[1:4]+param.est[9]*rep(udk1.pred,4)) /
  (1+exp(ddf[, 4:7]%*%param.est[1:4]+param.est[9]*rep(udk1.pred,4)) + exp(ddf[, 4:7]%*%param.est[5:8]+param.est[10]*rep(udk2.pred,4)))

pdj_unemployed <- exp(ddf[, 4:7]%*%param.est[5:8]+param.est[10]*rep(udk2.pred,4)) /
  (1+exp(ddf[, 4:7]%*%param.est[1:4]+param.est[9]*rep(udk1.pred,4)) + exp(ddf[, 4:7]%*%param.est[5:8]+param.est[10]*rep(udk2.pred,4)))

pdj_inact <-	1- pdj_employed - pdj_unemployed

# Census prediction
pred_censo <- Censo_2021T1$Conteo*cbind(pdj_employed, pdj_unemployed, pdj_inact)

# Sample data
pred_sample <- data.frame(aggregate(data.frame(lab1, lab2), by=list(PROV, SEXO1, age4), mean))
names(pred_sample) <- c('PROV', 'SEXO1', 'age4', 'lab1', 'lab2')
pred_sample$lab1 <- ndk*pred_sample$lab1
pred_sample$lab2 <- ndk*pred_sample$lab2

# Final prediction
predIN_lab1 <- (pred_censo[,1] + pred_sample$lab1) / Ndk
predIN_lab2 <- (pred_censo[,2] + pred_sample$lab2) / Ndk
predIN_lab3 <- 1 - predIN_lab1 - predIN_lab2

UNEMPRatef <- 100*(predIN_lab2 / (predIN_lab1 + predIN_lab2))


# Bootstrap resampling
B <- 500 # number of simulations

ydj_boot <- ydj_boot_age <- list()
pred_censo_real <- matrix(ncol=3, nrow=416)
realboot_lab1 <- realboot_lab2 <- realboot_lab3 <- predINboot_lab1 <- predINboot_lab2 <- predINboot_lab3 <-
  UNEMPRatef_real <- UNEMPRatef_boot <- matrix(ncol=416, nrow=B)

HajekUEMPRate <- Hajek1 <- Hajek2 <- Hajek3 <- matrix(ncol=416, nrow=B)

param.est_boot <- matrix(ncol=10, nrow=B)

for(b in 1:B){
  print(b)
  # Paso 1
  ## (a). Bootstrap population
  udk1.pred_boot <- rnorm(D, mean=0, 1)
  udk2.pred_boot <- rnorm(D, mean=0, 1)
  
  pdj_employed_boot <- exp(ddf[, 4:7]%*%param.est[1:4]+param.est[9]*rep(udk1.pred_boot,4)) /
                       (1 + exp(ddf[, 4:7]%*%param.est[1:4]+param.est[9]*rep(udk1.pred_boot,4)) + 
                        exp(ddf[, 4:7]%*%param.est[5:8]+param.est[10]*rep(udk2.pred_boot,4)))
  
  pdj_unemployed_boot <- exp(ddf[, 4:7]%*%param.est[5:8]+param.est[10]*rep(udk2.pred_boot,4)) / 
                         (1 + exp(ddf[,4:7]%*%param.est[1:4]+param.est[9]*rep(udk1.pred_boot,4)) +
                          exp(ddf[, 4:7]%*%param.est[5:8]+param.est[10]*rep(udk2.pred_boot,4)))
  
  pdj_inact_boot <- 1 - pdj_employed_boot - pdj_unemployed_boot
  
  # Simulated target variables for the sample
  # We generate the non-sample individuals for each cross dt (domain and age group).
  for(d in 1:D) {
    ydj_boot[[d]] <- ydj_boot_age[[d]] <- vector()
    # Determine the age group of the jth individual
    for(j in 1:nd[d]) {
      # age4-2
      if(tail(which(X[[d]][j,]==1),1)==2){
        ydj_boot[[d]] <- rbind(ydj_boot[[d]], t(rmultinom(1, size=1,
                                                          prob=c(pdj_employed_boot[d+104], pdj_unemployed_boot[d+104], pdj_inact_boot[d+104]))))
        ydj_boot_age[[d]] <- c(ydj_boot_age[[d]], 2)	
        # age4-3		
      } else if (tail(which(X[[d]][j,]==1),1)==3){
        ydj_boot[[d]] <- rbind(ydj_boot[[d]], t(rmultinom(1, size=1,
                                                          prob=c(pdj_employed_boot[d+208], pdj_unemployed_boot[d+208], pdj_inact_boot[d+208]))))
        ydj_boot_age[[d]] <- c(ydj_boot_age[[d]], 3)		
        # age4-4		
      } else if (tail(which(X[[d]][j,]==1),1)==4){
        ydj_boot[[d]] <- rbind(ydj_boot[[d]], t(rmultinom(1, size=1,
                                                          prob=c(pdj_employed_boot[d+312], pdj_unemployed_boot[d+312], pdj_inact_boot[d+312]))))
        ydj_boot_age[[d]] <- c(ydj_boot_age[[d]], 4)		
      } else{
        # age4-1
        ydj_boot[[d]] <- rbind(ydj_boot[[d]], t(rmultinom(1, size=1,
                                                          prob=c(pdj_employed_boot[d], pdj_unemployed_boot[d], pdj_inact_boot[d]))))
        ydj_boot_age[[d]] <- c(ydj_boot_age[[d]], 1)		
      }
    }
  }
  
  lab1boot <- unlist(sapply(ydj_boot, function(x) x[,1]))
  lab2boot <- unlist(sapply(ydj_boot, function(x) x[,2]))
  lab3boot <- unlist(sapply(ydj_boot, function(x) x[,3]))
  
  pred_sample <- data.frame(aggregate(data.frame(lab1boot, lab2boot), by=list(PROV, SEXO1, age4), mean))
  names(pred_sample) <- c('PROV', 'SEXO1', 'age4', 'lab1boot', 'lab2boot')
  pred_sample$lab1boot <- ndk*pred_sample$lab1boot
  pred_sample$lab2boot <- ndk*pred_sample$lab2boot
  
  # Results for the Hajek Direct estimator
  Hajek1[b, ] <- dir2(data=lab1boot, FACTOREL, domain=list(PROV,SEXO1, age4))[4]$mean
  
  Hajek2[b, ] <- dir2(data=lab2boot, FACTOREL, domain=list(PROV,SEXO1, age4))[4]$mean
  
  Hajek3[b, ] <- dir2(data=lab3boot, FACTOREL, domain=list(PROV,SEXO1, age4))[4]$mean
  
  HajekUEMPRate[b, ] <- 100 * Hajek2[b,] / (Hajek1[b,] + Hajek2[b,])
  
  # Bootstrap fitting
  mle.L_boot <- MLE_Laplace(x1, x2, ydj_boot, nud, seeds, seeds.u, max.iter=1000, epsilon.u=0.001, epsilon.par=0.0001)   
  
  param.est_boot[b,] <- t(as.matrix(mle.L_boot$par.est))
  udk.pred_boot <- mle.L_boot$udk.pred	
  
  udk1.pred_boot <- udk2.pred_boot <-rep(NA, D)
  for(d in 1:D) {
    udk1.pred_boot[d] <- udk.pred_boot[[d]][1]
    udk2.pred_boot[d] <- udk.pred_boot[[d]][2]
  }
  
  pdj_employed_boot2 <- exp(ddf[, 4:7]%*%param.est_boot[b,1:4]+param.est_boot[b,9]*rep(udk1.pred_boot,4)) /
                        (1+exp(ddf[, 4:7]%*%param.est_boot[b,1:4]+param.est_boot[b,9]*rep(udk1.pred_boot,4)) +
                        exp(ddf[, 4:7]%*%param.est_boot[b,5:8]+param.est_boot[b,10]*rep(udk2.pred_boot,4)))
  
  pdj_unemployed_boot2 <- exp(ddf[, 4:7]%*%param.est_boot[b,5:8]+param.est_boot[b,10]*rep(udk2.pred_boot,4)) /
                          (1+exp(ddf[, 4:7]%*%param.est_boot[b,1:4]+param.est_boot[b,9]*rep(udk1.pred_boot,4)) +
                          exp(ddf[, 4:7]%*%param.est_boot[b,5:8]+param.est_boot[b,10]*rep(udk2.pred_boot,4)))
  
  pdj_inact_boot2 <- 1- pdj_employed_boot2 - pdj_unemployed_boot2
  
  for (i in 1:416){
    pred_censo_real[i,] <- as.vector(rmultinom(1, size=Censo_2021T1$Conteo[i], 
                                               prob=c(pdj_employed_boot[i], pdj_unemployed_boot[i], pdj_inact_boot[i])))
  }
  
  realboot_lab1[b, ] <- (pred_censo_real[,1] + pred_sample$lab1boot) / Ndk
  realboot_lab2[b, ] <- (pred_censo_real[,2] + pred_sample$lab2boot) / Ndk
  realboot_lab3[b, ] <- 1 - realboot_lab1[b,]  - realboot_lab2[b,] 
  
  UNEMPRatef_real[b, ] <- 100*(realboot_lab2[b,] /(realboot_lab1[b,]  + realboot_lab2[b,] ))
  
  # Estimation
  pred_censo_boot <- Censo_2021T1$Conteo*cbind(pdj_employed_boot2, pdj_unemployed_boot2, pdj_inact_boot2)
  
  predINboot_lab1[b, ] <- (pred_censo_boot[,1] + pred_sample$lab1boot) / Ndk
  predINboot_lab2[b, ] <- (pred_censo_boot[,2] + pred_sample$lab2boot) / Ndk
  predINboot_lab3[b, ] <- 1 - predINboot_lab1[b,] - predINboot_lab2[b,]
  
  UNEMPRatef_boot[b, ] <- 100*predINboot_lab2[b,] / (predINboot_lab1[b,] + predINboot_lab2[b,])
  
}

# IC boot for the model parameters
print(apply(param.est_boot, 2, quantile, prob=c(0.025, 0.975)))

## Unemployment rates
# Results PI
rmse_UNEMPRatePI <- sqrt(colMeans((UNEMPRatef_real - UNEMPRatef_boot)^2))
rrmse_UNEMPRatePI <- rmse_UNEMPRatePI / UNEMPRatef

# Hajek (model-based approach)
rmse_UNEMPRateHajek <- sqrt(colMeans((UNEMPRatef_real - HajekUEMPRate)^2))
rrmse_UNEMPRateHajek <- rmse_UNEMPRateHajek / UNEMPRatef

save_results_UEMPRate <- data.frame(rmse_UNEMPRatePI, rrmse_UNEMPRatePI, rmse_UNEMPRateHajek, rrmse_UNEMPRateHajek)

write.table(save_results_UEMPRate, "mse_UNEMPRate_500.txt", sep=";", col.names = TRUE, row.names = FALSE)	

## Proportion of employed people
# Results PI
rmse_lab1PI <- sqrt(colMeans((realboot_lab1-predINboot_lab1)^2))
rrmse_lab1PI <- rmse_lab1PI / predIN_lab1

# Hajek (model-based approach)
rmse_lab1Hajek <- sqrt(colMeans((realboot_lab1-Hajek1)^2))
rrmse_lab1Hajek <- rmse_lab1Hajek / predIN_lab1

save_results_lab1 = data.frame(rmse_lab1PI, rrmse_lab1PI, rmse_lab1Hajek, rrmse_lab1Hajek )

write.table(save_results_lab1, "mse_lab1_500.txt", sep=";", col.names = TRUE, row.names = FALSE)

##  Proportion of unemployed people
# Results PI
rmse_lab2PI <- sqrt(colMeans((realboot_lab2-predINboot_lab2)^2))
rrmse_lab2PI <- rmse_lab2PI / predIN_lab2

# Hajek (model-based approach)
rmse_lab2Hajek <- sqrt(colMeans((realboot_lab2-Hajek2)^2))
rrmse_lab2Hajek <- rmse_lab2Hajek / predIN_lab2

save_results_lab2 <- data.frame(rmse_lab2PI, rrmse_lab2PI, rmse_lab2Hajek, rrmse_lab2Hajek)

write.table(save_results_lab2,"mse_lab2_500.txt", sep=";", col.names = TRUE, row.names = FALSE)

## Proportion of inactive people
# Results PI
rmse_lab3PI <- sqrt(colMeans((realboot_lab3-predINboot_lab3)^2))
rrmse_lab3PI <- rmse_lab3PI / predIN_lab3

# Hajek (model-based approach)
rmse_lab3Hajek <- sqrt(colMeans((realboot_lab3-Hajek3)^2))
rrmse_lab3Hajek <- rmse_lab3Hajek / predIN_lab3

save_results_lab3 <- data.frame(rmse_lab3PI, rrmse_lab3PI, rmse_lab3Hajek, rrmse_lab3Hajek)

write.table(save_results_lab3,"mse_lab3_500.txt", sep=";", col.names = TRUE, row.names = FALSE)
