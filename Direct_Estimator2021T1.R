



## This R code corresponds to the section "Application to real data" of the paper 
## Small area estimation of labour force indicators under unit-level multinomial mixed models.
## In this script we calculate the Hájek estimator, their variance and variation coefficient.


## AUTHORS: Bugallo M., Esteban M.D., Hozba T., Morales D., Perez A.



rm(list=ls())
library(car)

source('AuxFun.R') # Script of some auxiliary functions that must be in the same folder than this one.


########################################
######## PART 1 
########################################

# Read the csv file (must be in the same folder than this one)
EPA_2021T1 <- read.csv("EPA_2021T1.csv", sep=';', header=T)

# Delete individuals younger than 16
EPA_2021T1_new <- EPA_2021T1[which(EPA_2021T1$EDAD1>15),] 

attach(EPA_2021T1_new)
n <- nrow(EPA_2021T1_new)

# AOI: target variable
table(AOI)
lab <- recode(AOI, " c('3','4') = 'employed'; c('5','6') = 'unemployed'; c('7','8','9') = 'inactive' ", 
		as.factor=TRUE, levels = c('employed','unemployed', 'inactive'))

table(lab)
lab1 <- as.numeric(lab=='employed') # employed
lab2 <- as.numeric(lab=='unemployed') # unemployed
lab3 <- as.numeric(lab=='inactive') # inactive

age4 <- cut(EDAD1, breaks=c(0,45,55,64,Inf), include.lowest=TRUE, labels=c("1","2","3","4"))


########################################
######## PART 2 
########################################

## HAJEK DIRECT ESTIMATOR OF THE MEAN AND TOTAL #####
## GROUP BY PROVINCE, SEX AND AGE4

## Horvitz and Thompson for totals
dir <- data.frame(aggregate(FACTOREL*data.frame(1/FACTOREL, 1, lab1, lab2, lab3), by=list(PROV,SEXO1,age4), sum))
names(dir) <- c('Province', 'Sex', 'age4', 'nd', 'hatNd', 'lab1Total', 'lab2Total', 'lab3Total')

## Hajek for means
dir2.ds <- data.frame(Province=dir$Province, Sex=dir$Sex, age4=dir$age4, nd=dir$nd, hatNd=dir$hatNd)

# Estimates of means of employed people
dir2.ds$lab1mean <- dir$lab1Total/dir$hatNd

# Estimates of means of unemployed people
dir2.ds$lab2mean <- dir$lab2Total/dir$hatNd

# Estimates of means of inactive people
dir2.ds$lab3mean <- dir$lab3Total/dir$hatNd

## Variance of the Hajek direct estimator
dir2.ds$lab1meanvar <- dir2(data=lab1, FACTOREL, domain=list(PROV,SEXO1,age4))$var.mean
dir2.ds$lab2meanvar <- dir2(data=lab2, FACTOREL, domain=list(PROV,SEXO1,age4))$var.mean
dir2.ds$lab3meanvar <- dir2(data=lab3, FACTOREL, domain=list(PROV,SEXO1,age4))$var.mean

# CVs of the Hajek direct estimator
dir2.ds$lab1cv <- sqrt(dir2.ds$lab1meanvar)/abs(dir2.ds$lab1mean)
dir2.ds$lab1cv[dir2.ds$lab1mean==0] <- NaN

dir2.ds$lab2cv <- sqrt(dir2.ds$lab2meanvar)/abs(dir2.ds$lab2mean)
dir2.ds$lab2cv[dir2.ds$lab2mean==0] <- NaN

dir2.ds$lab3cv <- sqrt(dir2.ds$lab3meanvar)/abs(dir2.ds$lab3mean)
dir2.ds$lab3cv[dir2.ds$lab3mean==0] <- NaN

# Covariance between Hajek for means of employed and unemployed people
difference1 <- difference2 <- ww1 <- list() 
for(d in 1:nrow(dir2.ds)){
	condition <- paste(EPA_2021T1_new$PROV,EPA_2021T1_new$SEXO1,age4,sep="")==paste(dir2.ds$Province,dir2.ds$Sex,dir2.ds$age4,sep="")[d]
	difference1[[d]] <- lab1[condition]-dir2.ds$lab1mean[d] 
	difference2[[d]] <- lab2[condition]-dir2.ds$lab2mean[d] 
	ww1[[d]] <- FACTOREL[condition]*(FACTOREL[condition]-1)
}

ww1s1s2 <- mapply(ww1, mapply(difference1, difference2, FUN="*"), FUN="*")
dir2.ds$covarlab12 <- sapply(ww1s1s2, sum)/dir2.ds$hatNd^2

## Unemployment rate (in %)
dir2.ds$UNEMPRate <- 100 * dir2.ds$lab2mean/(dir2.ds$lab1mean+dir2.ds$lab2mean)

# Variance of the unemployment rate estimator
# (pp.55 eq 3.10 of A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL book available at www.springer.com)
s1.ds <- dir2.ds$lab2mean^2*dir2.ds$lab1meanvar/(dir2.ds$lab1mean+dir2.ds$lab2mean)^4 
s2.ds <- dir2.ds$lab1mean^2*dir2.ds$lab2meanvar/(dir2.ds$lab1mean+dir2.ds$lab2mean)^4 
s12.ds <- 2*dir2.ds$lab1mean*dir2.ds$lab2mean*dir2.ds$covarlab12/(dir2.ds$lab1mean+dir2.ds$lab2mean)^4 

dir2.ds$vrate <- 10^4 * (s1.ds+s2.ds-s12.ds)

# CVs of the unemployment rate estimator
dir2.ds$cvrate <- sqrt(dir2.ds$vrate)/abs(dir2.ds$UNEMPRate)
dir2.ds$cvrate[dir2.ds$UNEMPRate == 0 ] <- NaN

write.csv(dir2.ds, 'EPA_2021T1_Hajek.csv', row.names = FALSE)