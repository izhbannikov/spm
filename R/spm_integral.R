library(stats)


spm_integral <- function(dat,parameters) {
  
  res_prev <- NULL
  pars_prev <<- parameters
  iteration <- 0
  
  lower_bound <<- c(-1,0,1e-16,0,0,0, 0)
  upper_bound <<- c(0, Inf, 1e-6, Inf, Inf, Inf, 0.1) #c(0,100,1e-3,10,100,1, 1)
  
  
  maxlik <- function(dat, par) {
    a <- par[1] 
    f1 <- par[2]
    Q <- par[3]
    b <- par[4]
    f <- par[5]
    mu0 <- par[6]
    theta <- par[7]
    
    dims <- dim(dat)
    res <- .Call("complik", dat, dims[1], dims[2], a, f1, Q, b, f, mu0, theta)
    
    cat("L=",res,"\n")
    cat("Iter:", iteration, 
        "a=",par[1],
        "f1=",par[2],
        "Q=",par[3],
        "b=",par[4],
        "f=",par[5],
        "mu0=",par[6], 
        "theta=", par[7], "\n")
    
    iteration <<- iteration + 1
    
    # Adjusting lower and upper bounds:
    for(i in length(par)) {
      if(par[i] <= lower_bound[i]) {
        lower_bound[i] <<- lower_bound[i] - parameters[i]/3
        cat("Lower bound\n")
        cat(lower_bound)
      } else if(par[i] >= upper_bound[i]) {
        upper_bound[i] <<- upper_bound[i] + parameters[i]/3
        cat("Upper bound\n")
        cat(upper_bound)
      }
    }
    
    #if(is.nan(res)) {
    #  res <- maxlik(dat, pars_prev)
    ##  parameters <<- pars_prev
    #} else {
    #  res_prev <<- res
    #  pars_prev <<- par
    #}
    res
  }
  
  
  maxlik(as.matrix(dat), pars_prev)
  # Optimization:
  result <- optim(par = pars_prev, 
                  fn=maxlik, dat = as.matrix(dat), control = list(fnscale=-1, trace=T, maxit=10000, factr=1e-16, ndeps=c(1e-12,1e-12,1e-16,1e-12,1e-12,1e-12,1e-12)), 
                  method="L-BFGS-B", lower = lower_bound, upper = upper_bound)
  
  result
}

