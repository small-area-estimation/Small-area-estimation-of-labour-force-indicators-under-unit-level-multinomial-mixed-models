



## This R code corresponds to the section "Application to real data" of the paper 
## Small area estimation of labour force indicators under unit-level multinomial mixed models.
## In this script we implements Newton-Raphson and MLE-Laplace fit algorithms.

## AUTHORS: Bugallo M., Esteban M.D., Hozba T., Morales D., Perez A.




##############################
## Newton Raphson algorithm ##
##############################

NewtonRaphson_u <- function(param, X1=x1[[d]], X2=x2[[d]], y=ydj[[d]], nud.vec=nud[[d]], epsilon, maxiter, ud.seed)
{
  
  beta1 <- param[1:p1]
  beta2 <- param[(p1+1):p]
  phi.0 <- param[(p+1):(p+2)]
  udk.0 <- ud.seed    
  
  Sd <- matrix(0,2,1)
  Hd <- matrix(0,2,2)
  
  num.iter <- 0
  convergence <- FALSE
  
  while (convergence == FALSE & num.iter<=maxiter)
  {
    
    etad.1 <- X1 %*% beta1 + phi.0[1] * udk.0[1]
    etad.2 <- X2 %*% beta2 + phi.0[2] * udk.0[2]
    
    aux <- exp(etad.1) + exp(etad.2)
    
    pd.NR <- matrix(0,length(etad.1),2)
    
    pd.NR[,1] <- exp(etad.1)/(1+aux)
    pd.NR[,2] <- exp(etad.2)/(1+aux)
    
    
    for (k in 1:2){
      Sd[k,1] <- -phi.0[k] * nud.vec %*% pd.NR[,k] + phi.0[k] * sum(y[,k]) - udk.0[k]
      Hd[k,k] <- -1 - phi.0[k]^2 * nud.vec %*% (pd.NR[,k]*(1-pd.NR[,k]))
    }
    Hd[1,2] <- Hd[2,1] <- phi.0[1] * phi.0[2] * nud.vec %*% (pd.NR[,1]*pd.NR[,2])
    
    udk.1 <- udk.0 - solve(Hd) %*% Sd
    
    diff <- abs(udk.1-udk.0)
    
    if (any(diff>=epsilon)){
      udk.0 <- udk.1
      num.iter  <- num.iter+1
    }
    else{
      convergence <- TRUE
      num.iter  <- num.iter+1
    }  
  }
  return(list(u=udk.1, num.iter, convergence))
}     


##################################
## Laplace approximation of the ## 
##        log-likelihood        ##
##################################

log.lik.L <- function(param, x1=x1, x2=x2, ydj=ydj, nud = nud, udk.NR=udk.NR){
  
  beta1 <- param[1:p1]
  beta2 <- param[(p1+1):p]
  phi.0 <- param[(p+1):(p+2)]
  
  Hd <- matrix(0,2,2)
  
  sum.1 <- vector()
  sum.2 <- vector()
  term.3 <- vector()
  term.4 <- vector()
  
  
  for(d in 1:D){
    
    X1 <- x1[[d]]
    X2 <- x2[[d]]
    y <- ydj[[d]]
    nud.vec <- nud[[d]]
    udk.0 <- udk.NR[[d]]  
    
    #  calculating the matrix G_0d
    
    etad.1 <- X1 %*% beta1 + phi.0[1] * udk.0[1]
    etad.2 <- X2 %*% beta2 + phi.0[2] * udk.0[2]
    
    aux <- exp(etad.1) + exp(etad.2)
    
    pd.L <- matrix(0,length(etad.1),2)
    
    pd.L[,1] <- exp(etad.1)/(1+aux)
    pd.L[,2] <- exp(etad.2)/(1+aux)
    
    for (k in 1:2){
      Hd[k,k] <- -1 - phi.0[k]^2 * nud.vec %*% (pd.L[,k]*(1-pd.L[,k]))
    }
    Hd[1,2] <- Hd[2,1] <- phi.0[1] * phi.0[2] * nud.vec %*% (pd.L[,1]*pd.L[,2])
    
    Gd <- -Hd
    
    sum.1[d] <- nud.vec %*% log(1+aux)
    sum.2[d] <- y[,1] %*% etad.1 + y[,2] %*% etad.2
    
    term.3[d] <- (udk.0 %*% udk.0)/2 
    term.4[d] <- log(det(Gd))/2
    
  }
  
  return(- sum( -sum.1 + sum.2 - term.3 - term.4 ) )    
  
}


############################
##  MLE-Laplace algorithm ##
############################

MLE_Laplace <- function(x1, x2, ydj, nud, seeds.param, seeds.u, max.iter=1000, epsilon.u=0.001, epsilon.par=0.0001){
  
  param.0 <- seeds.param    # seeds for beta, phi
  
  udk.NR.0 <- udk.NR <- seeds.u      # seeds for u_dk
  
  num.iter <- 0
  convergence <- FALSE
  
  
  while (convergence == FALSE & num.iter<=max.iter)
  {
    
    #  Newton-Raphson for u_dk
    for (d in 1:D) {
      
      res <- NewtonRaphson_u(param=param.0, X1=x1[[d]], X2=x2[[d]], y=ydj[[d]], nud.vec=nud[[d]], epsilon=0.001, maxiter=100, ud.seed = udk.NR.0[[d]])
      
      udk.NR[[d]] <- as.vector(res[[1]])
      
    }
    
    res.opt <- nloptr(param.0, eval_f = log.lik.L, lb = c(rep(-Inf, p), rep(0, q-1)), ub = c(rep(Inf, p+q-1)), opts = opts,
                      x1=x1, x2=x2, ydj=ydj, nud = nud, udk.NR=udk.NR)
    
    param.1 <- res.opt$solution
    
    #### Here we check the difference between iterations
    
    u.NR.0 <- Reduce(udk.NR.0, f=rbind)
    u.NR <- Reduce(udk.NR, f=rbind)
    
    diff.u <- abs(u.NR - u.NR.0)
    diff.par <- abs(param.1 - param.0)
    
    if (any(diff.u>=epsilon.u) & any(diff.par>=epsilon.par)){
      udk.NR.0 <- udk.NR
      param.0 <- param.1
      num.iter  <- num.iter+1
    }
    else{
      convergence <- TRUE
      num.iter  <- num.iter+1
    }
    
  } 
  
  return(list(par.est = param.1, udk.pred = udk.NR, num.iter = num.iter, convergence = convergence))
  
}
