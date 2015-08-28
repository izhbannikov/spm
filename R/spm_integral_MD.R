library(stats)


spm_integral_MD <- function(dat,parameters) {
  
  #res_prev <- NULL
  #pars_prev <<- parameters
  #iteration <- 0
  #
  #lower_bound <<- c(-1,0,1e-15,0,0,0, 0)
  #upper_bound <<- c(0, Inf, 1e-13, Inf, Inf, Inf, 0.1) #c(0,100,1e-3,10,100,1, 1)
  
  maxlik <- function(dat, par) {
    a <- par$aH
    f1 <- par$f1H
    Q <- par$QH
    b <- par$bH
    f <- par$fH
    mu0 <- par$mu0
    theta <- par$theta
    k <- par$k
    #print(a)
    dims <- dim(dat)
    print(dims)
    #res <- .Call("complik", dat, dims[1], dims[2], a, f1, Q, b, f, mu0, theta)
    res <- .Call("complikMD", dat, dims[1], dims[2], a, f1, Q, b, f, mu0, theta, k)
    
    #cat("L=",res,"\n")
    #cat("Iter:", iteration, 
    #            "a=",par[1],
    #            "f1=",par[2],
    #            "Q=",par[3],
    #            "b=",par[4],
    #            "f=",par[5],
    #            "mu0=",par[6], 
    #            "theta=", par[7], "\n")
    #
    #iteration <<- iteration + 1
  
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

 
  maxlik(dat,parameters)
  # Optimization:
  #result <- optim(par = pars_prev, 
  #              fn=maxlik, dat = as.matrix(dat), control = list(fnscale=-1, trace=T, maxit=10000, factr=1e-16, ndeps=c(1e-8,1e-3,1e-12,1e-3,1e-3,1e-8,1e-8)), 
  #              method="L-BFGS-B", lower = lower_bound, upper = upper_bound)
  #
  #result
}


