library(stats)

results <- list(aH=NULL, f1H=NULL, QH=NULL, fH=NULL, bH=NULL, mu0H=NULL, thetaH=NULL)

setBoundaries <- function(k, params) {
  # This function sets lower and upper boundaries for optim.
  # - k - number of dimensions
  # - params - initial parameters, a vector
  #
  # Lower and upper boundaries:
  lower_bound <- c(); upper_bound <- c()
  # Setting boundaries for coefficients:
  # aH
  start=1; end=k^2
  lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){res=-1.5})))
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){res=-1e-12})))
  # f1H
  start=end+1; end=start+k-1
  lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] <= 0, 0.25*params[n], -0.25*params[n]) }))) 
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] >= 0, 0.25*params[n], -0.25*params[n]) })))
  # QH
  start=end+1; end=start+k^2-1
  lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n) {x=ifelse(params[n] >= 0, runif(1,0,1e-9), -1*runif(1,0,1e-7))})) )
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n) {x=ifelse(params[n] >= 0, runif(1,0,1e-7), -1*runif(1,0,1e-9))})))
  # fH
  start=end+1; end=start+k-1
  lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] <= 0, 0.25*params[n], -0.25*params[n]) })) )
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] >= 0, 0.25*params[n], -0.25*params[n]) })))
  # bH
  start=end+1; end=start+k-1
  lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] <= 0, 0.5*params[n], -0.5*params[n]) })) )
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] >= 0, 0.5*params[n], -0.5*params[n]) })))
  # mu0
  start=end+1; end=start
  lower_bound <- c(lower_bound, 1e-8 )
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] >= 0, 5*params[n], -5*params[n]) })))
  # theta
  start=end+1; end=start
  lower_bound <- c(lower_bound, 1e-6 )
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] >= 0, params[n], -1*params[n]) })))
  
  
  res=list(lower_bound=lower_bound, upper_bound=upper_bound)
}

spm_integral_MD <- function(dat,parameters, k, dd) {
  final_res <- list()
  iteration <- 0
  bounds <- setBoundaries(k, parameters)
  
  ndeps <- c(rep(1e-6,k^2), # a
              rep(1e-3,k), # f1
              rep(1e-16,k^2), # Q
              rep(1e-3,k), # f
              rep(1e-3,k), # b
              1e-6, # mu0
              1e-4) # theta
  
  
  maxlik <- function(dat, par) {
    stopflag <- F
    # Reading parameters:
    start=1; end=k^2
    a <- matrix(par[start:end],ncol=k, byrow=F); results$a <<- a
    start=end+1; end=start+k-1
    f1 <- matrix(par[start:end],ncol=k, byrow=F); results$f1 <<- a
    start=end+1; end=start+k^2-1
    Q <- matrix(par[start:end],ncol=k, byrow=F); results$Q <<- Q
    start=end+1; end=start+k-1
    b <- matrix(par[start:end],nrow=k); results$b <<- b
    start=end+1; end=start+k-1
    f <- matrix(par[start:end],ncol=k, byrow=F); results$f <<- f
    start=end+1; end=start
    mu0 <- par[start:end]; results$mu0 <<- mu0
    start=end+1; end=start
    theta <- par[start:end]; results$theta <<- theta
    # End reading parameters
    
    for(i in 1:length(results)) {
      if(length(intersect(results[[i]],c(bounds$lower_bound[i], bounds$upper_bound[i]))) >= 1) {
        print(results[[i]])
        stopflag <- T
      }
    }
    
    if(stopflag == F) {
      dims <- dim(dat)
      res <- .Call("complikMD", dat, dims[1], dims[2], a, f1, Q, b, f, mu0, theta, k)
    } else {
      cat("Optimization stopped. Parametes achieved lower or upper bound, you need more data to correctrly obtain optimal parameters.")
      return
    }
    
    iteration <<- iteration + 1
    
    cat("L = ",res,"\n")
    cat("Iteration: ", iteration,  "\nResults:\n") 
    print(results)
    
    res
  }

  # Optimization:
  optim_results <- NULL
  optim_results <- optim(par = parameters, 
                fn=maxlik, dat = as.matrix(dat), control = list(fnscale=-1, trace=T, factr=1e-16, ndeps=ndeps), 
                method="L-BFGS-B", lower = bounds$lower_bound, upper = bounds$upper_bound)
  
  final_res <<- list(results, optim_results)
  dd$results <<- final_res
  final_res
}


