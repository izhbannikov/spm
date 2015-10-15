library(stats)

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
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){res=0.2})))
  # f1H
  start=end+1; end=start+k-1
  lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){ 0 }))) 
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] > 0, 0.05*params[n], -0.05*params[n]) })))
  # QH
  start=end+1; end=start+k^2-1
  lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n) {x=ifelse(params[n] >= 0, runif(1,0,1e-12), -1*runif(1,0,1e-7))})) )
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n) {x=ifelse(params[n] > 0, runif(1,0,1e-7), -1*runif(1,0,1e-12))})))
  # fH
  start=end+1; end=start+k-1
  lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){ 0 })) )
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] > 0, 0.05*params[n], -0.05*params[n]) })))
  # bH
  start=end+1; end=start+k-1
  lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] <= 0, 0.5*params[n], -0.5*params[n]) })) )
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] > 0, 0.5*params[n], -0.5*params[n]) })))
  # mu0
  start=end+1; end=start
  lower_bound <- c(lower_bound, 1e-8 )
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] > 0, 5*params[n], -5*params[n]) })))
  # theta
  start=end+1; end=start
  lower_bound <- c(lower_bound, 1e-6 )
  upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] > 0, params[n], -1*params[n]) })))
  
  
  res=list(lower_bound=lower_bound, upper_bound=upper_bound)
}

#'Continuous multi-dimensional optimization
#'It is much slower that discrete but more precise and can handle time intervals with different lengths.
#'@param dat A data table.
#'@param parameters A starting pont (a vector).
#'@param k A number of dimensions.
#'@param verbose An indicator of verbosing output.
#'@param tol A tolerance threshold for matrix inversion.
#'@return A list of two elements: (1) parameters a, f1, Q, f, b, mu0, theta; (2) An output from "optim" function used 
#'for maximum likelihood estimation.
#'@examples
#'#'library(spm)
#'# Reading longitude data:
#'longdat <- read.csv(system.file("data","longdat.csv",package="spm"))
#'# Prepare data for optimization:
#'vitstat <- read.csv(system.file("data","vitstat.csv",package="spm"))
#'# Remove unneeded NAs:
#'longdat.nonan <- longdat[which(is.na(longdat$Age) == F),]
#'vitstat.nonan <- vitstat[which(is.na(vitstat$BirthCohort) == F),]
#'dat=prepare_data(longdat=longdat.nonan, vitstat=vitstat.nonan,interval=1, col.status="IsDead", col.id="ID", col.age="Age", col.age.next="AgeNext", col.age.event="LSmort", covariates=c("DBP"), verbose=T)
#'# Parameters estimation:
#'dat<-[,1:6]
#'pars=spm_integral_MD(dat, parameters=c(-0.05, 80, 2e-8, 80, 5, 2e-5, 0.08), k = 1)
#'pars
spm_integral_MD <- function(dat,parameters, k, verbose=F) {
  final_res <- list()
  # Current results:
  results <- list(a=NULL, f1=NULL, Q=NULL, f=NULL, b=NULL, mu0=NULL, theta=NULL)
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
    f1 <- matrix(par[start:end],ncol=k, byrow=F); results$f1 <<- f1
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
        cat("Parameter", names(results)[i], "achieved lower/upper bound. Process stopped.\n")
        cat(results[[i]],"\n")
        stopflag <- T
        break
      }
    }
    
    if(stopflag == F) {
      dims <- dim(dat)
      res <- .Call("complikMD", dat, dims[1], dims[2], a, f1, Q, b, f, mu0, theta, k)
      assign("results", results, envir=.GlobalEnv)
      iteration <<- iteration + 1
      if(verbose) {
        cat("L = ", res,"\n")
        cat("Iteration: ", iteration,  "\nResults:\n") 
        print(results)
      }
      
    } else {
      cat("Optimization stopped. Parametes achieved lower or upper bound.\nPerhaps you need more data or these returned parameters might be enough.\n")
      print("###########################################################")
      res <- NA
    }
    
    res
  }

  # Optimization:
  optim_results <- NA
  tryCatch(optim(par = parameters, 
                fn=maxlik, dat = as.matrix(dat), control = list(fnscale=-1, trace=T, factr=1e-16, ndeps=ndeps), 
                method="L-BFGS-B", lower = bounds$lower_bound, upper = bounds$upper_bound), 
           error=function(e) e, 
           finally=NA)
  final_res <<- list(results, optim_results)
  final_res
}


