library(stats)


spm_integral_MD <- function(dat,parameters) {
  
  res_prev <- NULL
  pars_prev <<- parameters[1:(length(parameters)-1)]
  iteration <- 0
  kk <- parameters[length(parameters)]
  
  setBoundaries <- function() {
    # Lower and upper boundaries:
    lower_bound <- c()
    upper_bound <- c()
    #
    start=1
    end=kk^2
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] <= 0, 2*pars_prev[n], -2*pars_prev[n]) })))
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] >= 0, 5*pars_prev[n], -5*pars_prev[n]) })))
    #
    start=end+1
    end=start+kk-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] <= 0, 5*pars_prev[n], -5*pars_prev[n]) })))
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] >= 0, 5*pars_prev[n], -5*pars_prev[n]) })))
    #
    start=end+1
    end=start+kk^2-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n) {x=ifelse(pars_prev[n] >= 0, 1e-12, -1e-7)})) )
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n) {x=ifelse(pars_prev[n] >= 0, 1e-6, -1e-12)})))
    #
    start=end+1
    end=start+kk-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] <= 0, 5*pars_prev[n], -5*pars_prev[n]) })) )
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] >= 0, 5*pars_prev[n], -5*pars_prev[n]) })))
    #
    start=end+1
    end=start+kk-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] <= 0, 5*pars_prev[n], -5*pars_prev[n]) })) )
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] >= 0, 5*pars_prev[n], -5*pars_prev[n]) })))
    # mu0
    start=end+1
    end=start
    lower_bound <- c(lower_bound, 1e-8 )
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] >= 0, 5*pars_prev[n], -5*pars_prev[n]) })))
    # theta
    start=end+1
    end=start
    lower_bound <- c(lower_bound, 1e-5 )
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){pars_prev[n] + ifelse(pars_prev[n] >= 0, pars_prev[n], -1*pars_prev[n]) })))
    
    res=list(lower_bound=lower_bound, upper_bound=upper_bound)
  }
  
  ndeps <- c(rep(1e-8,kk^2),
              rep(1e-4,kk),
              rep(1e-16,kk^2),
              rep(1e-4,kk),
              rep(1e-4,kk),
              1e-8,
              1e-6)
  
  
  maxlik <- function(dat, par) {
    
    # Adjusting lower and upper bounds:
    #for(i in length(par)) {
    #  if(par[i] <= lower_bound[i]) {
    #    lower_bound[i] <<- lower_bound[i] - parameters[i]/3
    #    cat("Lower bound\n")
    #    cat(lower_bound)
    #  } else if(par[i] >= upper_bound[i]) {
    #    upper_bound[i] <<- upper_bound[i] + parameters[i]/3
    #    cat("Upper bound\n")
    #    cat(upper_bound)
    #  }
    #}
    
    start=1
    end=kk^2
    a=matrix(par[start:end],ncol=kk, byrow=F)
    start=end+1
    end=start+kk-1
    f1 <- matrix(par[start:end],ncol=kk, byrow=F)
    start=end+1
    end=start+kk^2-1
    Q <- matrix(par[start:end],ncol=kk, byrow=F)
    start=end+1
    end=start+kk-1
    b <- matrix(par[start:end],nrow=kk)
    start=end+1
    end=start+kk-1
    f <- matrix(par[start:end],ncol=kk, byrow=F)
    start=end+1
    end=start
    mu0 <- par[start:end]
    start=end+1
    end=start
    theta <- par[start:end]
    
    dims <- dim(dat)
    res <- .Call("complikMD", dat, dims[1], dims[2], a, f1, Q, b, f, mu0, theta, kk)
    
    cat("L=",res,"\n")
    cat("Iter:", iteration, 
                "a=",a,
                "f1=",f1,
                "Q=",Q,
                "b=",b,
                "f=",f,
                "mu0=",mu0, 
                "theta=", theta, "\n")
    
    iteration <<- iteration + 1
  
    res
  }

  bounds <- setBoundaries()
  # Optimization:
  result <- optim(par = pars_prev, 
                fn=maxlik, dat = as.matrix(dat), control = list(fnscale=-1, trace=T, maxit=10000, factr=1e-16, ndeps=ndeps), 
                method="L-BFGS-B", lower = bounds$lower_bound, upper = bounds$upper_bound)
  
  result
}


