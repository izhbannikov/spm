library(stats)


spm_integral_MD <- function(dat,parameters) {
  
  res_prev <- NULL
  pars_prev <<- parameters[1:(length(parameters)-1)]
  iteration <- 0
  kk <- parameters[length(parameters)]
  #print(length(parameters))
  #print(parameters)
  lower_bound <<- c(rep(-1, kk^2),
                    rep(0,kk),
                    rep(1e-15, kk^2),
                    rep(0,kk),
                    rep(0,kk),
                    0, 
                    0)
  upper_bound <<- c(rep(0,kk^2),
                    rep(Inf,kk),
                    rep(1e-6,kk^2), 
                    rep(Inf,kk),
                    rep(Inf,kk),
                    1, 
                    0.1)
  
  ndeps <<- c(rep(1e-12,kk^2),
              rep(1e-12,kk),
              rep(1e-16,kk^2),
              rep(1e-12,kk),
              rep(1e-12,kk),
              1e-12,
              1e-12)
  
  #lower_bound <<- c(-1,0,1e-15,0,0,0, 0)
  #upper_bound <<- c(0, Inf, 1e-6, Inf, Inf, Inf, 0.1) #c(0,100,1e-3,10,100,1, 1)
  
  
  
  
  maxlik <- function(dat, par) {
    start=1
    end=kk^2
    a=matrix(par[start:end],ncol=kk, byrow=F)
    print(a)
    start=end+1
    end=start+kk-1
    f1 <- matrix(par[start:end],ncol=kk, byrow=F)
    print(f1)
    start=end+1
    end=start+kk^2-1
    Q <- matrix(par[start:end],ncol=kk, byrow=F)
    print(Q)
    start=end+1
    end=start+kk-1
    b <- matrix(par[start:end],nrow=kk)
    print(b)
    start=end+1
    end=start+kk-1
    f <- matrix(par[start:end],ncol=kk, byrow=F)
    print(f)
    start=end+1
    end=start
    mu0 <- par[start:end]
    print(mu0)
    start=end+1
    end=start
    theta <- par[start:end]
    print(theta)
    
    
    dims <- dim(dat)
    #print(dims)
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
    
    #if(is.nan(res)) {
    #  res <- maxlik(dat, pars_prev)
    ##  parameters <<- pars_prev
    #} else {
    #  res_prev <<- res
    #  pars_prev <<- par
    #}
    res
  }

  # Optimization:
  #print(pars_prev)
  #print(ndeps)
  result <- optim(par = pars_prev, 
                fn=maxlik, dat = as.matrix(dat), control = list(fnscale=-1, trace=T, maxit=10000, factr=1e-16, ndeps=ndeps), 
                method="L-BFGS-B", lower = lower_bound, upper = upper_bound)
  
  result
}


