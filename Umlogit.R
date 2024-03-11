



## This R code corresponds to the section "Application to real data" of the paper 
## Small area estimation of labour force indicators under unit-level multinomial mixed models.
## This is the main script that fits the Multinomial Mixed Model to the dataset EPA_2021T1.csv.


## AUTHORS: Bugallo M., Esteban M.D., Hozba T., Morales D., Perez A.



rm(list=ls())
source("function_MLE_Laplace.R") # Script that implements algorithms that must be in the same folder than this one.
source('AuxFun.R') # Script of some auxiliary functions that must be in the same folder than this one.

library(car)
library(nnet)
library(mefa)
library(lme4)
library(mvtnorm)
library(cubature)
library(nloptr)
library(fastDummies)
library(numDeriv)
library(KernSmooth)


## A unit-level multinomial logit mixed model by province:sex
# Fitting model

## EPA_2021T1
EPA_2021T1 <- read.csv("EPA_2021T1.csv", sep=';', header=T) 
EPA_2021T1_new <- EPA_2021T1[which(EPA_2021T1$EDAD1>15), ] 
attach(EPA_2021T1_new)

# Hajek results
EPA_2021T1_Hajek <- read.csv("EPA_2021T1_Hajek.csv", sep=',', header=T) # This file was created with Direct_Estimator2021T1.R script

# Target vector variable
lab <- recode(AOI, " c('3','4') = 'employed'; c('5','6') = 'unemployed'; c('7','8','9') = 'inactive' ", 
              as.factor=TRUE, levels = c('employed','unemployed', 'inactive'))
lab1 <- as.numeric(lab=='employed')
lab2 <- as.numeric(lab=='unemployed')
lab3 <- as.numeric(lab=='inactive')

# Auxiliary variables
age4 <- cut(EDAD1, breaks=c(0,45,55,64,Inf), include.lowest=TRUE, labels=c("1","2","3","4"))

# Census
Censo_2021T1 <- read.csv('Censo_2021T1.csv', sep=',', header=T)

n <- length(lab1)
D <- length(unique(PROV))*length(unique(SEXO1))

nd <- data.frame(aggregate(rep(1, n), by=list(PROV, SEXO1), sum))$x
ndk <- data.frame(aggregate(rep(1, n), by=list(PROV, SEXO1, age4), sum))$x

Ndk <- Censo_2021T1$Conteo + ndk


########################################
######## PART 1 
########################################

# Multinomial model without random effects

modmax <- multinom(relevel(lab, ref="inactive") ~ age4)
X <- model.matrix(modmax)

coeficientes <- summary(modmax)$coefficients
pl_employed <- X%*%coeficientes[1,]
pl_unemployed <- X%*%coeficientes[2,]

prob_employed <- exp(pl_employed)/(1+exp(pl_employed)+exp(pl_unemployed))
prob_unemployed <- exp(pl_unemployed)/(1+exp(pl_employed)+exp(pl_unemployed))

df <- data.frame(rbind(c(1,0,0,0), c(1,1,0,0), c(1,0,1,0), c(1,0,0,1)))
colnames(df) <- c('intercept', 'age42', 'age43', 'age44')
rownames(df) <- c('G1', 'G2', 'G3', 'G4')

pdj0_employed <- exp(as.matrix(df)%*%coeficientes[1,])/(1+exp(as.matrix(df)%*%coeficientes[1,]) +
                                                          exp(as.matrix(df)%*%coeficientes[2,]))

pdj0_unemployed <- exp(as.matrix(df)%*%coeficientes[2,])/(1+exp(as.matrix(df)%*%coeficientes[1,]) +
                                                            exp(as.matrix(df)%*%coeficientes[2,]))

pred0_censo <- Censo_2021T1$Conteo*cbind(rep(pdj0_employed, each=104), rep(pdj0_unemployed, each=104))

pred0_sample <- data.frame(aggregate(data.frame(lab1, lab2), by=list(PROV, SEXO1, age4), mean))
names(pred0_sample) <- c('PROV', 'SEXO1', 'age4', 'lab1', 'lab2')
pred0_sample$lab1 <- ndk*pred0_sample$lab1
pred0_sample$lab2 <- ndk*pred0_sample$lab2

# Final prediction
pred0_lab1 <- (pred0_censo[,1] + pred0_sample$lab1)/Ndk
pred0_lab2 <- (pred0_censo[,2] + pred0_sample$lab2)/Ndk
pred0_lab3 <- 1 - pred0_lab1 - pred0_lab2

UNEMPRatef0 <- 100*(pred0_lab2 / (pred0_lab1 + pred0_lab2))

resultados_fixed_effect <- data.frame('PROV'= rep(1:52, 4), 'SEXO1'=c(rep(1,52), rep(6, 52)), 
                                      pred0_lab1, pred0_lab2, pred0_lab3, UNEMPRatef0)


########################################
######## PART 2 
########################################
# Multinomial model with random effects by province:sex

set.seed(19032005)

# Domain indexes
dd <- rep(NA, n)
for (i in 1:n){
  if (SEXO1[i]==1){
    dd[i] <- PROV[i]
  } else
    dd[i] <- PROV[i] + 52
}

# age4 (for lab1 and lab2)
X <- model.matrix(modmax) 
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

# Target variable
ydj <- ydj2 <- list()

for(d in 1:D) {
  ydj[[d]] <- vector()
  for(j in 1:nd[d]) {
    ydj[[d]] <- rbind(ydj[[d]], t(c(lab1[matriz$dd==d][j], lab2[matriz$dd==d][j], lab3[matriz$dd==d][j])))
  }
  ydj2[[d]] <- ydj[[d]][,c(1,2)]
}

# Fitting the final model using the nloptr function
nudj <- 1
beta1.0 <- c(0.81, 0.36, -0.98, 4.44) # Obtained at code line 72 (multinomial model without random effects)
beta2.0 <- c(-0.78, -0.01, -1.40, -5.94) # Obtained at code line 72 (multinomial model without random effects)
beta.0 <- c(beta1.0, beta2.0)

p1 <- length(beta1.0)
p2 <- length(beta2.0)
p <- p1 + p2

phi1.0 <- 1
phi2.0 <- 1
phi.vect <- c(phi1.0,phi2.0)
phi <- as.matrix(bdiag(phi1.0, phi2.0))
q <- nrow(phi) + 1

param.0 <- seeds <- c(beta.0, phi1.0, phi2.0)

nud <- list()
for(d in 1:D) {
  nud[[d]] <- rep(1,nd[d])
}

a <- proc.time()
seeds.u <- split(matrix(0, D, 2), 1:D)  
opts <- list( "algorithm" = "NLOPT_LN_NEWUOA_BOUND", maxeval=1000, print_level=3)    
mle.L <- MLE_Laplace(x1, x2, ydj, nud, seeds, seeds.u, max.iter=1000, epsilon.u=0.001, epsilon.par=0.0001) 
b <- proc.time()
b-a	

param.est <- as.matrix(mle.L$par.est)
num.iter <- mle.L$num.iter
convergence <- mle.L$convergence
udk.pred <- mle.L$udk.pred

udk1.pred <- sapply(udk.pred, function(x) x[1])
udk2.pred <- sapply(udk.pred, function(x) x[2])

# Asympothic Standar Error 
log.lik.L_exp <- function(param) { log.lik.L(param, x1=x1, x2=x2, ydj=ydj, nud=nud, udk.NR=udk.pred) }

hess_parametros <- solve(hessian(log.lik.L_exp, param.est))
se <- sqrt(diag(hess_parametros))

# p value
pvalue <- function(param, sd) {
  z <- abs(param)/sd
  p <- pnorm(z, lower.tail=F)
  return( 2*p )
}
round(pvalue(param.est, se), 2)

# IC
round(data.frame('LB'=param.est-qnorm(0.975)*se, 'UB'=param.est+qnorm(0.975)*se), 4)


########################################
######## PART 3 
########################################

### Plug-in probabilities

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
pdj_inact <- 1 - pdj_employed - pdj_unemployed

UNEMPRate <- 100 * pdj_unemployed / (pdj_employed + pdj_unemployed)


# Some plots included in the paper

par(mar=c(4, 4.3, 4, 4), xpd=T)
plot(EPA_2021T1_Hajek$UNEMPRate, type='l', col='gray40', main='Unemployment rates', ylab='Percentage', xaxt = "n", xlab='Age group', lwd=1.5, 
     cex.main=1.8, yaxt = "n", cex.lab=1.7, ylim=c(0,70)); 
axis(1, at=c(0, 52, 104, 156, 208, 260, 312, 364, 416), labels=c('', 1, '', 2, '', 3, '', 4, ''), cex.axis=1.7); 
axis(2, cex.axis=1.7); lines(UNEMPRate, type='l', col='red')	
legend('topleft', col=c('gray40','red'), c('Hajek', 'Plug-in'), bty = "n", cex=2, pt.cex = 2, lwd=2)


par(mar=c(4, 4.3, 4, 4), xpd=T)
plot(pdj_employed[EPA_2021T1_Hajek$Sex==1], EPA_2021T1_Hajek$lab1mean[EPA_2021T1_Hajek$Sex==1], pch=19, col='darkblue', main= 'Employed',
     ylab='Hajek', xlab='Plug-in', xlim=c(0,0.9), ylim=c(0,0.9), cex.main=1.8, cex.lab=1.8, cex.axis=1.8)
points(pdj_employed[EPA_2021T1_Hajek$Sex==6], EPA_2021T1_Hajek$lab1mean[EPA_2021T1_Hajek$Sex==6], pch=3, col='darkblue')	
new <- data.frame(x = seq(0, 0.9, by=0.01))
lines(new$x, new$x, col='red', lwd=2, lty=2)
fit0 <- locpoly(pdj_employed, EPA_2021T1_Hajek$lab1mean, bandwidth = 0.05, degree=1) 	
lines(fit0, lwd=3)
legend('topleft', col='darkblue', c('Men', 'Women'), bty = "n", cex=2, pt.cex = 1.5, pch=c(19,3))

par(mar=c(4, 4.3, 4, 4), xpd=T)
plot(pdj_unemployed[EPA_2021T1_Hajek$Sex==1], EPA_2021T1_Hajek$lab2mean[EPA_2021T1_Hajek$Sex==1], pch=19, col='darkblue', main= 'Unemployed',
     ylab='Hajek', xlab='Plug-in', xlim=c(0,0.3), ylim=c(0,0.3), cex.main=1.8, cex.lab=1.8, cex.axis=1.8)
points(pdj_unemployed[EPA_2021T1_Hajek$Sex==6], EPA_2021T1_Hajek$lab2mean[EPA_2021T1_Hajek$Sex==6],	pch=3, col='darkblue')	
new <- data.frame(x = seq(0, 0.3, by=0.01))
lines(new$x, new$x, col='red', lwd=2, lty=2)
fit0 <- locpoly(pdj_unemployed, EPA_2021T1_Hajek$lab2mean, bandwidth = 0.05, degree=1) 	
lines(fit0, lwd=3)

par(mar=c(4, 4.3, 4, 4), xpd=T)
plot(pdj_inact[EPA_2021T1_Hajek$Sex==1], EPA_2021T1_Hajek$lab3mean[EPA_2021T1_Hajek$Sex==1], pch=19, col='darkblue', main= 'Inactive',
     ylab='Hajek', xlab='Plug-in', xlim=c(0,1), ylim=c(0,1), cex.main=1.8, cex.lab=1.8, cex.axis=1.8)
points(pdj_inact[EPA_2021T1_Hajek$Sex==6], EPA_2021T1_Hajek$lab3mean[EPA_2021T1_Hajek$Sex==6], pch=3, col='darkblue')	
new <- data.frame(x = seq(0, 1, by=0.01))
lines(new$x, new$x, col='red', lwd=2, lty=2)
fit0 <- locpoly(pdj_inact, EPA_2021T1_Hajek$lab3mean, bandwidth = 0.06, degree=1) 	
lines(fit0, lwd=3)

par(mar=c(4, 4.3, 4, 4), xpd=T)
plot(UNEMPRate[EPA_2021T1_Hajek$Sex==1], EPA_2021T1_Hajek$UNEMPRate[EPA_2021T1_Hajek$Sex==1], pch=19, col='darkblue',
     main= 'Unemployment rates', ylab='Hajek (%)', xlab='Plug-in (%)', xlim=c(0,70), ylim=c(0,70), cex.main=1.8, cex.lab=1.8, cex.axis=1.8)
points(UNEMPRate[EPA_2021T1_Hajek$Sex==6], EPA_2021T1_Hajek$UNEMPRate[EPA_2021T1_Hajek$Sex==6], pch=3, col='darkblue')	
new <- data.frame(x = seq(0, 70, by=0.01))
lines(new$x, new$x, col='red', lwd=2, lty=2)
fit0 <- locpoly(UNEMPRate, EPA_2021T1_Hajek$UNEMPRate, bandwidth = 1.3, degree=1)
lines(fit0, lwd=3)


########################################
######## PART 4 
########################################

# Model validation

pdk1 <- data.frame(aggregate(lab1, by=list(PROV, SEXO1, age4), sum))$x/ndk
pdk2 <- data.frame(aggregate(lab2, by=list(PROV, SEXO1, age4), sum))$x/ndk
pdk3 <- data.frame(aggregate(lab3, by=list(PROV, SEXO1, age4), sum))$x/ndk

# Standarized residuals
ASRlab1<- (pdk1 - pdj_employed) / sd(pdk1 - pdj_employed)
ASRlab2<- (pdk2 - pdj_unemployed) / sd(pdk2 - pdj_unemployed)

# Boxplots
par(mar=c(4, 4.3, 4, 4), xpd=T)
boxplot(ASRlab1~EPA_2021T1_Hajek$Province, main='Province', ylab='', xlab='', cex.main=1.8, cex.lab=1.8, cex.axis=1.8)
EPA_2021T1_Hajek$Sex2=EPA_2021T1_Hajek$Sex
EPA_2021T1_Hajek$Sex2[EPA_2021T1_Hajek$Sex2==6]=2
boxplot(ASRlab1~EPA_2021T1_Hajek$Sex2, main='Sex', ylab='', xlab='', cex.main=1.8, cex.lab=1.8, cex.axis=1.8)
boxplot(ASRlab1~EPA_2021T1_Hajek$age4, main='Age group', ylab='', xlab='', cex.main=1.8, cex.lab=1.8, cex.axis=1.8)

boxplot(ASRlab2~EPA_2021T1_Hajek$Province, main='Province', ylab='', xlab='', cex.main=1.8, cex.lab=1.8, cex.axis=1.8)
boxplot(ASRlab2~EPA_2021T1_Hajek$Sex2, main='Sex', ylab='', xlab='', cex.main=1.8, cex.lab=1.8, cex.axis=1.8)
boxplot(ASRlab2~EPA_2021T1_Hajek$age4, main='Age group', ylab='', xlab='', cex.main=1.8, cex.lab=1.8, cex.axis=1.8)

# Other plots
par(mar=c(4, 4.3, 4, 4), xpd=T)	 
plot(ASRlab1[sort(ndk, index.return=T)$ix], type='l', ylab=' ', xlab='Subdomain sample size', cex.main=1.8, cex.lab=1.7,
     cex.axis=1.7, main='Employed', xaxt = 'n')
new <- data.frame(x = seq(1, 416))
lines(new$x, rep(0,416), col='red', lwd=2, lty=2)		
axis(1, at=round(seq(1,416,length=6),0), labels=floor(sort(ndk, index.return=T)$x)[round(seq(1,416,length=6),0)], cex.axis=1.7)
axis(2, cex.axis=1.7)

par(mar=c(4, 4.3, 4, 4), xpd=T)	 	
plot(ASRlab2[sort(ndk, index.return=T)$ix], type='l', ylab=' ',
     xlab='Subdomain sample size', cex.main=1.8, cex.lab=1.7, cex.axis=1.7, main='Unemployed', xaxt = 'n')
new <- data.frame(x = seq(1, 416))
lines(new$x, rep(0,416), col='red', lwd=2, lty=2)		
axis(1, at=round(seq(1,416,length=6),0), labels=floor(sort(ndk, index.return=T)$x)[round(seq(1,416,length=6),0)], cex.axis=1.7)
axis(2, cex.axis=1.7)


########################################
######## PART 7 
########################################

# Final prediction

# Census 
pred_censo <- Censo_2021T1$Conteo * cbind(pdj_employed, pdj_unemployed, pdj_inact)

pred_lab1 <- (pred_censo[,1] + pred0_sample$lab1)/Ndk
pred_lab2 <- (pred_censo[,2] + pred0_sample$lab2)/Ndk
pred_lab3 <- 1 - pred_lab1 - pred_lab2

UNEMPRatef <- 100*(pred_lab2 / (pred_lab1 + pred_lab2))

resultados_random_effect <- data.frame('PROV'= rep(1:52, 4), 'SEXO1'=c(rep(1,52), rep(6, 52)), pred_lab1, pred_lab2, pred_lab3, UNEMPRatef)	


# Saving results
# 1.
UNEMPRate_11 <- UNEMPRatef[EPA_2021T1_Hajek$age4==1 & EPA_2021T1_Hajek$Sex==1][1:50]
UNEMPRate_16 <- UNEMPRatef[EPA_2021T1_Hajek$age4==1 & EPA_2021T1_Hajek$Sex==6][1:50]

UNEMPRate_21 <- UNEMPRatef[EPA_2021T1_Hajek$age4==2 & EPA_2021T1_Hajek$Sex==1][1:50]
UNEMPRate_26 <- UNEMPRatef[EPA_2021T1_Hajek$age4==2 & EPA_2021T1_Hajek$Sex==6][1:50]

UNEMPRate_31 <- UNEMPRatef[EPA_2021T1_Hajek$age4==3 & EPA_2021T1_Hajek$Sex==1][1:50]
UNEMPRate_36 <- UNEMPRatef[EPA_2021T1_Hajek$age4==3 & EPA_2021T1_Hajek$Sex==6][1:50]

UNEMPRate_41 <- UNEMPRatef[EPA_2021T1_Hajek$age4==4 & EPA_2021T1_Hajek$Sex==1][1:50]
UNEMPRate_46 <- UNEMPRatef[EPA_2021T1_Hajek$age4==4 & EPA_2021T1_Hajek$Sex==6][1:50]

UNEMPL_RATE <- data.frame('prov'=EPA_2021T1_Hajek$Province[1:50], UNEMPRate_11, UNEMPRate_16, UNEMPRate_21, UNEMPRate_26,
                          UNEMPRate_31, UNEMPRate_36, UNEMPRate_41, UNEMPRate_46)
write.table(UNEMPL_RATE, "UNEMP.RATE.txt", sep=";", col.names = TRUE, row.names = FALSE)


# 2.
inac_11 <- pred_lab3[EPA_2021T1_Hajek$age4==1 & EPA_2021T1_Hajek$Sex==1][1:50]
inac_16 <- pred_lab3[EPA_2021T1_Hajek$age4==1 & EPA_2021T1_Hajek$Sex==6][1:50]

inac_21 <- pred_lab3[EPA_2021T1_Hajek$age4==2 & EPA_2021T1_Hajek$Sex==1][1:50]
inac_26 <- pred_lab3[EPA_2021T1_Hajek$age4==2 & EPA_2021T1_Hajek$Sex==6][1:50]

inac_31 <- pred_lab3[EPA_2021T1_Hajek$age4==3 & EPA_2021T1_Hajek$Sex==1][1:50]
inac_36 <- pred_lab3[EPA_2021T1_Hajek$age4==3 & EPA_2021T1_Hajek$Sex==6][1:50]

inac_41 <- pred_lab3[EPA_2021T1_Hajek$age4==4 & EPA_2021T1_Hajek$Sex==1][1:50]
inac_46 <- pred_lab3[EPA_2021T1_Hajek$age4==4 & EPA_2021T1_Hajek$Sex==6][1:50]	

INACTIVE <- data.frame('prov'=EPA_2021T1_Hajek$Province[1:50], inac_11, inac_16, inac_21, inac_26, inac_31, inac_36, inac_41, inac_46)
write.table(INACTIVE, "INACTIVE.txt", sep=";", col.names = TRUE, row.names = FALSE)


# 3.
ocup_11 <- pred_lab1[EPA_2021T1_Hajek$age4==1 & EPA_2021T1_Hajek$Sex==1][1:50]
ocup_16 <- pred_lab2[EPA_2021T1_Hajek$age4==1 & EPA_2021T1_Hajek$Sex==6][1:50]

ocup_21 <- pred_lab1[EPA_2021T1_Hajek$age4==2 & EPA_2021T1_Hajek$Sex==1][1:50]
ocup_26 <- pred_lab1[EPA_2021T1_Hajek$age4==2 & EPA_2021T1_Hajek$Sex==6][1:50]

ocup_31 <- pred_lab1[EPA_2021T1_Hajek$age4==3 & EPA_2021T1_Hajek$Sex==1][1:50]
ocup_36 <- pred_lab1[EPA_2021T1_Hajek$age4==3 & EPA_2021T1_Hajek$Sex==6][1:50]

ocup_41 <- pred_lab1[EPA_2021T1_Hajek$age4==4 & EPA_2021T1_Hajek$Sex==1][1:50]
ocup_46 <- pred_lab1[EPA_2021T1_Hajek$age4==4 & EPA_2021T1_Hajek$Sex==6][1:50]	

EMPLOYED <- data.frame('prov'=EPA_2021T1_Hajek$Province[1:50], ocup_11, ocup_16, ocup_21, ocup_26, ocup_31, ocup_36, ocup_41, ocup_46)
write.table(EMPLOYED, "EMPLOYED.txt", sep=";", col.names = TRUE, row.names = FALSE)


# 4.
parado_11 <- pred_lab2[EPA_2021T1_Hajek$age4==1 & EPA_2021T1_Hajek$Sex==1][1:50]
parado_16 <- pred_lab2[EPA_2021T1_Hajek$age4==1 & EPA_2021T1_Hajek$Sex==6][1:50]

parado_21 <- pred_lab2[EPA_2021T1_Hajek$age4==2 & EPA_2021T1_Hajek$Sex==1][1:50]
parado_26 <- pred_lab2[EPA_2021T1_Hajek$age4==2 & EPA_2021T1_Hajek$Sex==6][1:50]

parado_31 <- pred_lab2[EPA_2021T1_Hajek$age4==3 & EPA_2021T1_Hajek$Sex==1][1:50]
parado_36 <- pred_lab2[EPA_2021T1_Hajek$age4==3 & EPA_2021T1_Hajek$Sex==6][1:50]

parado_41 <- pred_lab2[EPA_2021T1_Hajek$age4==4 & EPA_2021T1_Hajek$Sex==1][1:50]
parado_46 <- pred_lab2[EPA_2021T1_Hajek$age4==4 & EPA_2021T1_Hajek$Sex==6][1:50]	


UNEMPLOYED <- data.frame('prov'=EPA_2021T1_Hajek$Province[1:50], parado_11, parado_16, parado_21, parado_26, parado_31, parado_36,
                         parado_41, parado_46)
write.table(UNEMPLOYED, "UNEMPLOYED.txt", sep=";", col.names = TRUE, row.names = FALSE)	



