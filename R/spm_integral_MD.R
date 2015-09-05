library(stats)


spm_integral_MD <- function(dat,parameters, k, dd) {
  final_res <<- list()
  setBoundaries <- function(kk) {
    # Lower and upper boundaries:
    lower_bound <- c()
    upper_bound <- c()
    # aH
    start=1
    end=kk^2
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){res=-1.5})))
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){res=-1e-12})))
    # f1H
    start=end+1
    end=start+kk-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] <= 0, 0.25*pars_prev[n], -0.25*pars_prev[n]) })))
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] >= 0, 0.25*pars_prev[n], -0.25*pars_prev[n]) })))
    # QH
    start=end+1
    end=start+kk^2-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n) {x=ifelse(pars_prev[n] >= 0, runif(1,0,1e-9), -1*runif(1,0,1e-7))})) )
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n) {x=ifelse(pars_prev[n] >= 0, runif(1,0,1e-7), -1*runif(1,0,1e-9))})))
    # fH
    start=end+1
    end=start+kk-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] <= 0, 0.25*pars_prev[n], -0.25*pars_prev[n]) })) )
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] >= 0, 0.25*pars_prev[n], -0.25*pars_prev[n]) })))
    # bH
    start=end+1
    end=start+kk-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] <= 0, 0.5*pars_prev[n], -0.5*pars_prev[n]) })) )
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] >= 0, 0.5*pars_prev[n], -0.5*pars_prev[n]) })))
    # mu0
    start=end+1
    end=start
    lower_bound <- c(lower_bound, 1e-8 )
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] >= 0, 5*pars_prev[n], -5*pars_prev[n]) })))
    # theta
    start=end+1
    end=start
    lower_bound <- c(lower_bound, 1e-6 )
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] >= 0, pars_prev[n], -1*pars_prev[n]) })))
    
    res=list(lower_bound=lower_bound, upper_bound=upper_bound)
  }
  
  results <- list(aH=NULL, f1H=NULL, QH=NULL, fH=NULL, bH=NULL, mu0H=NULL, thetaH=NULL)
  pars_prev <<- parameters
  iteration <- 0
  bounds <<- setBoundaries(k)
  
  ndeps <- c(rep(1e-6,k^2),
              rep(1e-3,k),
              rep(1e-16,k^2),
              rep(1e-3,k),
              rep(1e-3,k),
              1e-6,
              1e-4)
  
  
  maxlik <- function(dat, par) {
    stopflag <- F
    # Reading parameters:
    start=1
    end=k^2
    a <- matrix(par[start:end],ncol=k, byrow=F) 
    results$aH <<- a
    start=end+1
    end=start+k-1
    f1 <- matrix(par[start:end],ncol=k, byrow=F)
    results$f1H <<- a
    start=end+1
    end=start+k^2-1
    Q <- matrix(par[start:end],ncol=k, byrow=F)
    results$QH <<- Q
    start=end+1
    end=start+k-1
    b <- matrix(par[start:end],nrow=k)
    results$bH <<- b
    start=end+1
    end=start+k-1
    f <- matrix(par[start:end],ncol=k, byrow=F)
    results$fH <<- f
    start=end+1
    end=start
    mu0 <- par[start:end]
    results$mu0 <<- mu0
    start=end+1
    end=start
    theta <- par[start:end]
    results$thetaH <<- theta
    # End reading parameters
    
    opt_pars <- list(a=a, f1=f1, Q=Q, f=f, b=b, mu0=mu0, theta=theta)
    for(i in 1:length(opt_pars)) {
      if(length(intersect(opt_pars[[i]],c(bounds$lower_bound[i], bounds$upper_bound[i]))) >= 1) {
        print(opt_pars[[i]])
        stopflag <- T
      }
    }
    
    if(stopflag == F) {
      dims <- dim(dat)
      res <- .Call("complikMD", dat, dims[1], dims[2], a, f1, Q, b, f, mu0, theta, k)
    } else {
      cat("Optimization stopped. Parametes achieved lower or upper bound, you need more data to correctrly obtain optimal parameters.")
      assign(results, envir=dd)
      stop()
    }
    
    iteration <<- iteration + 1
    
    cat("L = ",res,"\n")
    cat("Iteration: ", iteration,  "\nResults:\n") 
    print(results)
    
    res
  }

  
  
  # Optimization:
  optim_results <- optim(par = pars_prev, 
                fn=maxlik, dat = as.matrix(dat), control = list(fnscale=-1, trace=T, maxit=10000, factr=1e-16, ndeps=ndeps), 
                method="L-BFGS-B", lower = bounds$lower_bound, upper = bounds$upper_bound)
  
  final_res <<- list(results, optim_results)
  dd$results <<- final_res
  final_res
}


