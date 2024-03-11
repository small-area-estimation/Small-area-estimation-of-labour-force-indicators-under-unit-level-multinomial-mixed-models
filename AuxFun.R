



## This R code corresponds to the section "Application to real data" of the paper 
## Small area estimation of labour force indicators under unit-level multinomial mixed models.
## This script containts some auxiliary functions.

## AUTHORS: Bugallo M., Esteban M.D., Hozba T., Morales D., Perez A.




############################
## Hajek direct estimator ##
############################

dir2 <- function(data, w, domain, Nd) { 
  if(is.vector(data)){
    last <- length(domain) + 1
    Nd.hat <- aggregate(w, by=domain, sum)[,last]
    nd <- aggregate(rep(1, length(data)), by=domain, sum)[,last] 
    Sum <- aggregate(w*data, by=domain, sum)
    mean <- Sum[,last]/Nd.hat
    dom <- as.numeric(Reduce("paste0", domain)) 
    if(length(domain)==1){
      domain.unique <- sort(unique(dom)) }
    else{
      domain.unique <- as.numeric(Reduce("paste0", Sum[,1:length(domain)]))
    }
    difference <- list() 
    for(d in 1:length(mean)){
      condition <- dom==domain.unique[d]
      difference[[d]] <- w[condition]*(w[condition]-1)*(data[condition]-mean[d])^2 }
    var.mean <- unlist(lapply(difference, sum))/Nd.hat^2 
    if(missing(Nd)){
      return(data.frame(Sum[,-last], mean, var.mean, Nd.hat, nd)) }
    else{
      tot <- mean*Nd
      var.tot <- var.mean*Nd^2
      return(data.frame(Sum[,-last], tot, var.tot, mean, var.mean, Nd.hat, Nd, nd))
    } }
  else{
    warning("Only a numeric or integer vector must be called as data", call. = FALSE)
  } 
}


#############################
##  Model-based predictors ##
##  Sections 4.3 and 4.4   ##
#############################

# Plug-in predictor

IN <- function(ydj2, dk.pred, param.est){
	
  phiud <- sapply(udk.pred, FUN="%*%", as.matrix(bdiag(param.est[9], param.est[10])))
  
  etad.11 <- param.est[1] + phiud[1,]
  etad.21 <- param.est[5] + phiud[2,]
  etad.12 <- param.est[1] + param.est[2] + phiud[1,]
  etad.22 <- param.est[5] + param.est[5] + phiud[2,]
  
  etad.13 <- param.est[1] + param.est[3] + phiud[1,]
  etad.23 <- param.est[5] + param.est[7] + phiud[2,]
  etad.14 <- param.est[1] + param.est[4] + phiud[1,]
  etad.24 <- param.est[5] + param.est[8]+ phiud[2,]
  
  qd.1.11 <- exp(etad.11) / (1+exp(etad.11)+exp(etad.21))
  qd.2.11 <- exp(etad.21) / (1+exp(etad.11)+exp(etad.21))
  qd.1.22 <- exp(etad.12) / (1+exp(etad.12)+exp(etad.22))
  qd.2.22 <- exp(etad.22) / (1+exp(etad.12)+exp(etad.22))
  
  qd.1.33 <- exp(etad.13) / (1+exp(etad.13)+exp(etad.23))
  qd.2.33 <- exp(etad.23) / (1+exp(etad.13)+exp(etad.23))
  qd.1.44 <- exp(etad.14) / (1+exp(etad.14)+exp(etad.24))
  qd.2.44 <- exp(etad.24) / (1+exp(etad.14)+exp(etad.24))
  
  y1mean.in <-(1/Nd)*(sapply(ydj2, colSums)[1,] + Ndr.11*qd.1.11 + Ndr.22*qd.1.22 + Ndr.33*qd.1.33 + Ndr.11*qd.1.44)
  y2mean.in <-(1/Nd)*(sapply(ydj2, colSums)[2,] + Ndr.11*qd.2.11 + Ndr.22*qd.2.22 + Ndr.33*qd.2.33 + Ndr.11*qd.2.44)
  
  Rd.in <- y2mean.in/(y1mean.in  + y2mean.in)

  return(list(y1mean.in, y2mean.in, Rd.in))	
}


# Non-NAN mean and auxiliary functions

mean2 <- function (x) sum(x*is.finite(x), na.rm=TRUE)/length(is.finite(x))

exp_1 <- function(x, suma){ return(exp(x+1400)) }
exp_2 <- function(x, suma){ return(exp(x+2400)) }
exp_3 <- function(x, suma){ return(exp(x+2900)) }
exp_4 <- function(x, suma){ return(exp(x+3400)) }
exp_5 <- function(x, suma){ return(exp(x+3900)) }

f_inf <- function(z){return(z[!is.infinite(z)])}
not_inf <- function (y) { lapply(y, FUN=f_inf) }
	
external <- function (x) {
	D.d.S1 = vector('list', d)
	means <- sapply(x , FUN=mean)
	
	a_1 = (1:d)[means > (-2500) & means < (-1300)]
	D.d.S1[a_1] = lapply(x[a_1], FUN=exp_1)
	
	a_2 = (1:d)[means > (-3000) & means < (-2500)]
	D.d.S1[a_2] = lapply(x[a_2], FUN=exp_2)
	
	a_3 = (1:d)[means > (-3500) & means < (-3000)]
	D.d.S1[a_3] = lapply(x[a_3], FUN=exp_3)
	
	a_4 = (1:d)[means > (-4000) & means < (-3500)]
	D.d.S1[a_4] = lapply(x[a_4], FUN=exp_4)
	
	a_5 = (1:d)[means > (-4500) & means < (-4000)]
	D.d.S1[a_5] = lapply(x[a_4], FUN=exp_5)
	
	a_6 = (1:104)[-c(a_1, a_2, a_3, a_4, a_5)]
	D.d.S1[a_6] = lapply(x[a_6], FUN=exp)
	
	D.d.S1 = not_inf(D.d.S1)
	
	return(D.d.S1)
	}
	

# EBP predictor

EBP <- function(x1, x2, ydj2, param.est, S1){	
	# Paso 3.9i
	udk.s <- list()
	for(d in 1:D){
    	udk.s[[d]] <- matrix(0,nrow=2*S1,ncol=2)
    	for(s in 1:S1){
     	 udk.s[[d]][s,] <- rnorm(2)
      	udk.s[[d]][s+S1,] <- -1*udk.s[[d]][s,]
		}  
	}
  
	# Paso 3.9ii
	phiud.s <- lapply(udk.s, FUN="%*%", as.matrix(bdiag(param.est[9], param.est[10])))
	
	etad.s.11 <- etad.s.21 <- etad.s.12 <- etad.s.22 <- etad.s.13 <- etad.s.23 <-etad.s.14 <- etad.s.24 <- list()
	
	etad.s.11 <- lapply(lapply(phiud.s, FUN=function(x) x[,1]), param.est[1], FUN="+")
	etad.s.21 <- lapply(lapply(phiud.s, FUN=function(x) x[,2]), param.est[5], FUN="+")
	etad.s.12 <- lapply(lapply(phiud.s, FUN=function(x) x[,1]), param.est[1] + param.est[2], FUN="+")
  	etad.s.22 <- lapply(lapply(phiud.s, FUN=function(x) x[,2]), param.est[5] + param.est[6], FUN="+")
  	
  	etad.s.13 <- lapply(lapply(phiud.s, FUN=function(x) x[,1]), param.est[1] + param.est[3], FUN="+")
  	etad.s.23 <- lapply(lapply(phiud.s, FUN=function(x) x[,2]), param.est[5] + param.est[7], FUN="+")
  	etad.s.14 <- lapply(lapply(phiud.s, FUN=function(x) x[,1]), param.est[1] + param.est[4], FUN="+")
  	etad.s.24 <- lapply(lapply(phiud.s, FUN=function(x) x[,2]), param.est[5] + param.est[8], FUN="+")
	
	
	den.11 <- Map(lapply(etad.s.11, FUN=exp), lapply(lapply(etad.s.21, FUN=exp), 1, FUN="+"), f="+")	
	
 	qd.s.1.11 <- Map(lapply(etad.s.11, FUN=exp), den.11, f="/") 
  	qd.s.2.11 <- Map(lapply(etad.s.21, FUN=exp), den.11, f="/")
  
    den.22 <- Map(lapply(etad.s.12, FUN=exp), lapply(lapply(etad.s.22, FUN=exp), 1, FUN="+"), f="+")
    qd.s.1.22 <- Map(lapply(etad.s.12, FUN=exp), den.22, f="/") 
    qd.s.2.22 <- Map(lapply(etad.s.22, FUN=exp), den.22, f="/")
   
    den.33 <- Map(lapply(etad.s.13, FUN=exp), lapply(lapply(etad.s.23, FUN=exp), 1, FUN="+"), f="+")
    qd.s.1.33 <- Map(lapply(etad.s.13, FUN=exp), den.33, f="/")
    qd.s.2.33 <- Map(lapply(etad.s.23, FUN=exp), den.33, f="/")
    
    den.44 <- Map(lapply(etad.s.14, FUN=exp), lapply(lapply(etad.s.24, FUN=exp), 1, FUN="+"), f="+")
    qd.s.1.44 <- Map(lapply(etad.s.14, FUN=exp), den.44, f="/")
    qd.s.2.44 <- Map(lapply(etad.s.24, FUN=exp), den.44, f="/")
    
    # Paso 3.9iii
  	sum1D.d <- lapply(lapply(lapply(Map(lapply(ydj2, FUN=colSums), lapply(udk.s, FUN=t), f="*"), FUN=t), FUN="%*%",
  	as.matrix(bdiag(param.est[9], param.est[10]))), FUN=rowSums)
  	
  	aux1 <- aux2 <- list() 
    x1.beta <- lapply(x1, FUN="%*%", param.est[1:4])
    phiud1.s <- lapply(udk.s, function(x) x[,1])
    x2.beta <- lapply(x2, FUN="%*%", param.est[5:8])
    phiud2.s <- lapply(udk.s, function(x) x[,2])
    
    
    for(d in 1:D){
    	cat("d=", d,"\n")
    	aux1[[d]] <- matrix(0, nrow=2*S1, ncol=nd[[d]])
    	aux2[[d]] <- matrix(0, nrow=2*S1, ncol=nd[[d]])
    	for(j in 1:nd[[d]]){
        	for(s in (1:(2*S1))){
        		aux1[[d]][s,j] <- x1.beta[[d]][j] + phiud1.s[[d]][s]
        		aux2[[d]][s,j] <- x2.beta[[d]][j] + phiud2.s[[d]][s]
      	}
    	}
  	}
    sum2D.d <- lapply(lapply(lapply(Map(lapply(aux1, FUN=exp), 
    	lapply(aux2, FUN=exp), f="+"), FUN="+", 1), FUN=log), FUN=rowSums)
    	
    	aux_D.d.S1 <- Map(sum1D.d, sum2D.d, f="-")
    		
    D.d.S1 <- external(aux_D.d.S1)
 	
  	D.d <- sapply(D.d.S1, FUN=mean2) # toma valores iguales a 0
  	
  	# Paso 3.9iv y 3.9v
  	
  	qd.s.1.11.ebp <- mapply(lapply(lapply(Map(qd.s.1.11, D.d.S1, f="*"), FUN=sum), 2*S1, FUN="/"), D.d, FUN="/")
  	qd.s.2.11.ebp <- mapply(lapply(lapply(Map(qd.s.2.11, D.d.S1, f="*"), FUN=sum), 2*S1, FUN="/"), D.d, FUN="/")
  	qd.s.1.22.ebp <- mapply(lapply(lapply(Map(qd.s.1.22, D.d.S1, f="*"), FUN=sum), 2*S1, FUN="/"), D.d, FUN="/")
  	qd.s.2.22.ebp <- mapply(lapply(lapply(Map(qd.s.2.22, D.d.S1, f="*"), FUN=sum), 2*S1, FUN="/"), D.d, FUN="/")
  	
  	qd.s.1.33.ebp <- mapply(lapply(lapply(Map(qd.s.1.33, D.d.S1, f="*"), FUN=sum), 2*S1, FUN="/"), D.d, FUN="/")
  	qd.s.2.33.ebp <- mapply(lapply(lapply(Map(qd.s.2.33, D.d.S1, f="*"), FUN=sum), 2*S1, FUN="/"), D.d, FUN="/")
  	qd.s.1.44.ebp <- mapply(lapply(lapply(Map(qd.s.1.44, D.d.S1, f="*"), FUN=sum), 2*S1, FUN="/"), D.d, FUN="/")
  	qd.s.2.44.ebp <- mapply(lapply(lapply(Map(qd.s.2.44, D.d.S1, f="*"), FUN=sum), 2*S1, FUN="/"), D.d, FUN="/")
  	
  	 # Paso 3.9vi
	y1mean.ebp <-(1/Nd)*(sapply(ydj2, colSums)[1,] + Ndr.11*qd.s.1.11.ebp + Ndr.22*qd.s.1.22.ebp + 
		Ndr.33*qd.s.1.33.ebp +Ndr.44*qd.s.1.44.ebp)
		
  	y2mean.ebp <-(1/Nd)*(sapply(ydj2, colSums)[2,] + Ndr.11*qd.s.2.11.ebp + Ndr.22*qd.s.2.22.ebp + 
  		Ndr.33*qd.s.2.33.ebp + Ndr.44*qd.s.2.44.ebp)
  	
  	Rd.in.ebp <- y2mean.ebp/(y1mean.ebp  + y2mean.ebp) 
  
  	# Paso 3.9vii  	
  	sum.yd.11 <- sum.yd.22 <- sum.yd.33 <- sum.yd.44 <- list()
  	sum.yd.11.s <- sum.yd.22.s <- sum.yd.33.s <- sum.yd.44.s <- list()
  	
  	# Para que todo este bien definido
  	qd.s.1.11.ebp[qd.s.1.11.ebp==0]=min(qd.s.1.11.ebp[qd.s.1.11.ebp!=0 & is.finite(qd.s.1.11.ebp) & is.na(qd.s.1.11.ebp)==FALSE])
  	qd.s.2.11.ebp[qd.s.2.11.ebp==0]=min(qd.s.2.11.ebp[qd.s.2.11.ebp!=0 & is.finite(qd.s.2.11.ebp) & is.na(qd.s.2.11.ebp)==FALSE])
  	qd.s.1.22.ebp[qd.s.1.22.ebp==0]=min(qd.s.1.22.ebp[qd.s.1.22.ebp!=0 & is.finite(qd.s.1.22.ebp) & is.na(qd.s.1.22.ebp)==FALSE])
  	qd.s.2.22.ebp[qd.s.2.22.ebp==0]=min(qd.s.2.22.ebp[qd.s.2.22.ebp!=0 & is.finite(qd.s.2.22.ebp) & is.na(qd.s.2.22.ebp)==FALSE])
  	
  	qd.s.1.33.ebp[qd.s.1.33.ebp==0]=min(qd.s.1.33.ebp[qd.s.1.33.ebp!=0 & is.finite(qd.s.1.33.ebp) & is.na(qd.s.1.33.ebp)==FALSE])
  	qd.s.2.33.ebp[qd.s.2.33.ebp==0]=min(qd.s.2.33.ebp[qd.s.2.33.ebp!=0 & is.finite(qd.s.2.33.ebp) & is.na(qd.s.2.33.ebp)==FALSE])
  	qd.s.1.44.ebp[qd.s.1.44.ebp==0]=min(qd.s.1.44.ebp[qd.s.1.44.ebp!=0 & is.finite(qd.s.1.44.ebp) & is.na(qd.s.1.44.ebp)==FALSE])
  	qd.s.2.44.ebp[qd.s.2.44.ebp==0]=min(qd.s.2.44.ebp[qd.s.2.44.ebp!=0 & is.finite(qd.s.2.44.ebp) & is.na(qd.s.2.44.ebp)==FALSE])
  	
  	qd.s.1.11.ebp[is.infinite(qd.s.1.11.ebp)]=max(qd.s.1.11.ebp[qd.s.1.11.ebp!=0 & is.finite(qd.s.1.11.ebp) & is.na(qd.s.1.11.ebp)==FALSE])
  	qd.s.2.11.ebp[is.infinite(qd.s.2.11.ebp)]=max(qd.s.2.11.ebp[qd.s.2.11.ebp!=0 & is.finite(qd.s.2.11.ebp) & is.na(qd.s.2.11.ebp)==FALSE])
  	qd.s.1.22.ebp[is.infinite(qd.s.1.22.ebp)]=max(qd.s.1.22.ebp[qd.s.1.22.ebp!=0 & is.finite(qd.s.1.22.ebp) & is.na(qd.s.1.22.ebp)==FALSE])
  	qd.s.2.22.ebp[is.infinite(qd.s.2.22.ebp)]=max(qd.s.2.22.ebp[qd.s.2.22.ebp!=0 & is.finite(qd.s.2.22.ebp) & is.na(qd.s.2.22.ebp)==FALSE])
  	
  	qd.s.1.33.ebp[is.infinite(qd.s.1.33.ebp)]=max(qd.s.1.33.ebp[qd.s.1.33.ebp!=0 & is.finite(qd.s.1.33.ebp) & is.na(qd.s.1.33.ebp)==FALSE])
  	qd.s.2.33.ebp[is.infinite(qd.s.2.33.ebp)]=max(qd.s.2.33.ebp[qd.s.2.33.ebp!=0 & is.finite(qd.s.2.33.ebp) & is.na(qd.s.2.33.ebp)==FALSE])
  	qd.s.1.44.ebp[is.infinite(qd.s.1.44.ebp)]=max(qd.s.1.44.ebp[qd.s.1.44.ebp!=0 & is.finite(qd.s.1.44.ebp) & is.na(qd.s.1.44.ebp)==FALSE])
  	qd.s.2.44.ebp[is.infinite(qd.s.2.44.ebp)]=max(qd.s.2.44.ebp[qd.s.2.44.ebp!=0 & is.finite(qd.s.2.44.ebp) & is.na(qd.s.2.44.ebp)==FALSE])
  	
  	qd.s.1.11.ebp[is.na(qd.s.1.11.ebp)]=median(qd.s.1.11.ebp[qd.s.1.11.ebp!=0 & is.finite(qd.s.1.11.ebp) & is.na(qd.s.1.11.ebp)==FALSE])
  	qd.s.2.11.ebp[is.na(qd.s.2.11.ebp)]=median(qd.s.2.11.ebp[qd.s.2.11.ebp!=0 & is.finite(qd.s.2.11.ebp) & is.na(qd.s.2.11.ebp)==FALSE])
  	qd.s.1.22.ebp[is.na(qd.s.1.22.ebp)]=median(qd.s.1.22.ebp[qd.s.1.22.ebp!=0 & is.finite(qd.s.1.22.ebp) & is.na(qd.s.1.22.ebp)==FALSE])
  	qd.s.2.22.ebp[is.na(qd.s.2.22.ebp)]=median(qd.s.2.22.ebp[qd.s.2.22.ebp!=0 & is.finite(qd.s.2.22.ebp) & is.na(qd.s.2.22.ebp)==FALSE])
  	
  	qd.s.1.33.ebp[is.na(qd.s.1.33.ebp)]=median(qd.s.1.33.ebp[qd.s.1.33.ebp!=0 & is.finite(qd.s.1.33.ebp) & is.na(qd.s.1.33.ebp)==FALSE])
  	qd.s.2.33.ebp[is.na(qd.s.2.33.ebp)]=median(qd.s.2.33.ebp[qd.s.2.33.ebp!=0 & is.finite(qd.s.2.33.ebp) & is.na(qd.s.2.33.ebp)==FALSE])
  	qd.s.1.44.ebp[is.na(qd.s.1.44.ebp)]=median(qd.s.1.44.ebp[qd.s.1.44.ebp!=0 & is.finite(qd.s.1.44.ebp) & is.na(qd.s.1.44.ebp)==FALSE])
  	qd.s.2.44.ebp[is.na(qd.s.2.44.ebp)]=median(qd.s.2.44.ebp[qd.s.2.44.ebp!=0 & is.finite(qd.s.2.44.ebp) & is.na(qd.s.2.44.ebp)==FALSE])
  	
  	for(s in 1:(2*S1)){
  	    cat("s=", s,"\n")
    		yd.11 <- yd.22 <- yd.33 <- yd.44 <- list()
    		yd.11.s <- yd.22.s <- yd.33.s <- yd.44.s <- list()
    	
    		for(d in 1:D){
     		yd.11[[d]] <- rmultinom(Nd.11[d]-nd.11[d], size=1, prob=c(qd.s.1.11.ebp[d], qd.s.2.11.ebp[d]))
      		yd.22[[d]] <- rmultinom(Nd.22[d]-nd.22[d], size=1, prob=c(qd.s.1.22.ebp[d], qd.s.2.22.ebp[d]))
     	 	yd.33[[d]] <- rmultinom(Nd.33[d]-nd.33[d], size=1, prob=c(qd.s.1.33.ebp[d], qd.s.2.33.ebp[d]))
      		yd.44[[d]] <- rmultinom(Nd.44[d]-nd.44[d], size=1, prob=c(qd.s.1.44.ebp[d], qd.s.2.44.ebp[d]))
      		
      		yd.11.s[[d]] <- matrix(ydj2[[d]][Map(lapply(x1, FUN=function(x) x[,2]==0), 
      			lapply(x2, FUN=function(x) x[,2]==0), f="&")[[d]],], ncol=2)
      		yd.22.s[[d]] <- matrix(ydj2[[d]][Map(lapply(x1, FUN=function(x) x[,2]==0), 
      			lapply(x2, FUN=function(x) x[,2]==0), f="&")[[d]],], ncol=2)
      		yd.33.s[[d]] <- matrix(ydj2[[d]][Map(lapply(x1, FUN=function(x) x[,2]==0), 
      			lapply(x2, FUN=function(x) x[,2]==0), f="&")[[d]],], ncol=2)
      		yd.44.s[[d]] <- matrix(ydj2[[d]][Map(lapply(x1, FUN=function(x) x[,2]==0), 
      			lapply(x2, FUN=function(x) x[,2]==0), f="&")[[d]],], ncol=2)
      			
      		}

	sum.yd.11[[s]] <- sapply(yd.11, FUN=rowSums)
    sum.yd.22[[s]] <- sapply(yd.22, FUN=rowSums)
    sum.yd.33[[s]] <- sapply(yd.33, FUN=rowSums)
    sum.yd.44[[s]] <- sapply(yd.44, FUN=rowSums)
    
    sum.yd.11.s[[s]] <- sapply(yd.11.s, FUN=colSums)
    sum.yd.22.s[[s]] <- sapply(yd.22.s, FUN=colSums)
    sum.yd.33.s[[s]] <- sapply(yd.33.s, FUN=colSums)
    sum.yd.44.s[[s]] <- sapply(yd.44.s, FUN=colSums)
    }
    
  	sum.j <- Map( lapply(Map(Map(Map(sum.yd.11, sum.yd.22, f="+"), sum.yd.33, f="+"), sum.yd.44, f="+"), 
  		function(x) x[-3,]),  Map(Map(Map(sum.yd.11.s, sum.yd.22.s, f="+"), sum.yd.33.s, f="+"), 
  			sum.yd.44.s, f="+"), f="+")
  	cociente <- Map( lapply(sum.j, FUN=function(x) x[2,]), lapply(sum.j, FUN=colSums), f="/")
  	Rd.ebp <- Reduce(cociente, f="+")/(2*S1)
  	
  	return(list(y1mean.ebp, y2mean.ebp, Rd.in.ebp, Rd.ebp))	
 } 	