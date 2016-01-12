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
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] >= 0, 2*params[n], -2*params[n]) })))
  # f1H
  start=end+1; end=start+k-1
  lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){ 0 }))) 
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] >= 0, 2*params[n], -0.5*params[n]) })))
  # QH
  start=end+1; end=start+k^2-1
  lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n) {params[n] - ifelse(params[n] >= 0, 10*params[n], -10*params[n] )})) )
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n) {params[n] + ifelse(params[n] > 0, 10*params[n], -10*params[n] )})))
  # fH
  start=end+1; end=start+k-1
  lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){ 0 })) )
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] > 0, 2*params[n], -2*params[n]) })))
  # bH
  start=end+1; end=start+k-1
  lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] <= 0, 2*params[n], -2*params[n]) })) )
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] > 0, 2*params[n], -2*params[n]) })))
  # mu0
  start=end+1; end=start
  lower_bound <- c( lower_bound, 0 )
  upper_bound <- c( upper_bound, 1 )
  # theta
  start=end+1; end=start
  lower_bound <- c( lower_bound, 1e-6 )
  upper_bound <- c( upper_bound, 1 )
  
  
  res=list(lower_bound=lower_bound, upper_bound=upper_bound)
}

#'Continuous multi-dimensional optimization
#'@param dat A data table.
#'@param a A starting value of the rate of adaptive response to any deviation of Y from f1(t).
#'@param f1 A starting value of the average age trajectories of the variables which process is forced to follow. 
#'@param Q Starting values of the quadratic hazard term.
#'@param f A starting value of the "optimal" value of variable which corresponds to the minimum of hazard rate at a respective time.
#'@param b A starting value of a diffusion coefficient representing a strength of the random disturbance from Wiener Process.
#'@param mu0 A starting value of the baseline hazard.
#'@param theta A starting value of the parameter theta (axe displacement of Gompertz function).
#'@param k A number of dimensions.
#'@param verbose An indicator of verbosing output.
#'@param tol A tolerance threshold for matrix inversion.
#'@return A set of estimated parameters a, f1, Q, f, b, mu0, theta.
#'@details \code{spm_continuous} runs much slower that discrete but more precise and can handle time intervals with different lengths.
#'@examples
#'library(spm)
#'# Reading the data:
#'dd <- prepare_data(x=read.csv(system.file("data","longdat.csv",package="spm")), y=read.csv(system.file("data","vitstat.csv",package="spm")))
#'data <- dd[[1]][,2:6] # We have to remove subject ID from the simulated data.
#'#Parameters estimation:
#'pars <- spm_continuous(dat=data,a=-0.05, f1=80, Q=2e-8, f=80, b=5, mu0=2e-5, theta=0.08, k = 1)
#'pars
spm_continuous <- function(dat, 
                           a=-0.05, 
                           f1=80, 
                           Q=2e-8,
                           f=81,
                           b=5,
                           mu0=2e-5,
                           theta=0.08,
                           k=1, 
                           stopifbound=FALSE, 
                           algorithm="NLOPT_LN_NELDERMEAD",
                           verbose=FALSE) {
  
  avail_algorithms <- c("NLOPT_LN_NEWUOA",
                        "NLOPT_LN_NEWUOA_BOUND",
                        "NLOPT_LN_NELDERMEAD")
  
  if(!(algorithm %in% avail_algorithms)) {
    stop(cat("Provided algorithm", algorithm, "not in the list of available optimization methods."))
  }
  
  dat <<- dat
  final_res <- list()
  
  if(mu0 < 0) {mu0 <- 0}
  
  parameters <- c(a, f1, Q, f, b, mu0, theta)
  # Current results:
  results <<- list(a=NULL, f1=NULL, Q=NULL, f=NULL, b=NULL, mu0=NULL, theta=NULL)
  results_tmp <<- list(a=NULL, f1=NULL, Q=NULL, f=NULL, b=NULL, mu0=NULL, theta=NULL)
  iteration <- 0
  bounds <- setBoundaries(k, parameters)
  
  ndeps <- c(rep(1e-6,k^2), # a
              rep(1e-3,k), # f1
              rep(1e-16,k^2), # Q
              rep(1e-3,k), # f
              rep(1e-3,k), # b
              1e-6, # mu0
              1e-4) # theta
  
  # Reading parameters:
  start=1; end=k^2
  a <- matrix(parameters[start:end],ncol=k, byrow=F)
  results$a <<- a
  start=end+1; end=start+k-1
  f1 <- matrix(parameters[start:end],ncol=k, byrow=F)
  results$f1 <<- f1
  start=end+1; end=start+k^2-1
  Q <- matrix(parameters[start:end],ncol=k, byrow=F)
  results$Q <<- Q
  start=end+1; end=start+k-1
  f <- matrix(parameters[start:end],ncol=k, byrow=F)
  results$f <<- f
  start=end+1; end=start+k-1
  b <- matrix(parameters[start:end],nrow=k)
  results$b <<- b
  start=end+1; end=start
  mu0 <- parameters[start:end]
  results$mu0 <<- mu0
  start=end+1; end=start
  theta <- parameters[start:end]
  results$theta <<- theta
  # End of reading parameters
  
  #maxlik <- function(dat, par) {
  maxlik <- function(par) {
    dat <<- dat
    stopflag <- F
    # Reading parameters:
    start=1; end=k^2
    a <- matrix(par[start:end],ncol=k, byrow=F)
    results_tmp$a <<- a
    start=end+1; end=start+k-1
    f1 <- matrix(par[start:end],ncol=k, byrow=F)
    results_tmp$f1 <<- f1
    start=end+1; end=start+k^2-1
    Q <- matrix(par[start:end],ncol=k, byrow=F)
    results_tmp$Q <<- Q
    start=end+1; end=start+k-1
    f <- matrix(par[start:end],ncol=k, byrow=F)
    results_tmp$f <<- f
    start=end+1; end=start+k-1
    b <- matrix(par[start:end],nrow=k)
    results_tmp$b <<- b
    start=end+1; end=start
    mu0 <- par[start:end]
    results_tmp$mu0 <<- mu0
    start=end+1; end=start
    theta <- par[start:end]
    results_tmp$theta <<- theta
    # End reading parameters
    
    if(stopifbound) {
      for(i in 1:length(results_tmp)) {
        if(length(intersect(results_tmp[[i]],c(bounds$lower_bound[i], bounds$upper_bound[i]))) >= 1) {
          cat("Parameter", names(results)[i], "achieved lower/upper bound. Process stopped.\n")
          cat(results_tmp[[i]],"\n")
          stopflag <- T
          break
        }
      }
    }
    
    if(stopflag == F) {
      dims <- dim(dat)
      res <<- .Call("complikMD", dat, dims[1], dims[2], a, f1, Q, b, f, mu0, theta, k)
      assign("results", results_tmp, envir=.GlobalEnv)
      iteration <<- iteration + 1
      if(verbose) {
        cat("L = ", res,"\n")
        cat("Iteration: ", iteration,  "\nResults:\n") 
        print(results_tmp)
      }
      
    } else {
      cat("Optimization stopped. Parametes achieved lower or upper bound.\nPerhaps you need more data or these returned parameters might be enough.\n")
      print("###########################################################")
      res <<- get("results",envir=.GlobalEnv)
    }
    
    res <- -1*res
    res
  }

  # Optimization:
  if(verbose) {
    cat("Lower bound:\n")
    print(bounds$lower_bound)
    cat("Upper bound:\n")
    print(bounds$upper_bound)
  }
  tryCatch(nloptr(x0 = parameters, 
                 eval_f = maxlik, opts = list("algorithm"=algorithm, 
                                              "xtol_rel"=1.0e-8),
                 lb = bounds$lower_bound, ub = bounds$upper_bound),  
           error=function(e) {if(verbose  == TRUE) {print(e)}}, 
           finally=NA)
  
  
  
  res <- get("results",envir=.GlobalEnv)
  invisible(res)
}


