#'Discrete multi-dimensional optimization
#'@param dat A data table.
#'@param k A number of dimensions.
#'@param theta_range A range of theta parameter (axe displacement of Gompertz function), default: from 0.001 to 0.09 with step of 0.001.
#'@param tol A tolerance threshold for matrix inversion (NULL by default).
#'@return A list of two elements: (1) estimated parameters u, R, b, Sigma, Q, mu0, theta and
#'(2) estimated parameters a, f1, Q, f, b, mu0, theta. Note: b and mu0 from first list are different 
#'from b and mu0 from the second list.
#'@details This function is way much faster that continuous \code{spm_continuous_MD(...)} (but less precise) and used mainly in 
#'estimation a starting point for the \code{spm_continuous_MD(...)}.
#'@examples
#'library(spm)
#'# Reading longitudinal data
#'longdat <- read.csv(system.file("data","longdat.csv",package="spm"))
#'# Prepare data for optimization
#'vitstat <- read.csv(system.file("data","vitstat.csv",package="spm"))
#'data <- prepare_data(longdat=longdat, vitstat=vitstat,interval=1, col.status="IsDead", col.id="ID", col.age="Age", col.age.event="LSmort", covariates=c("DBP"), verbose=T)
#'# Parameters estimation
#'pars <- spm_discrete(data[[2]], k=1, theta_range=seq(0.001,0.09,by=0.001), tol=NULL)
#'pars
spm_discrete2 <- function(dat,k=1, tol=NULL, parameters = c(mu0=5e-4, Q=6e-9, b=-1.15e-6, theta=0.01)) {
  # Call function that computes Likelihood:
  parameters <- parameters
  
  maxlik <- function(dat, par) {
    dims <- dim(dat)
    # Reading parameters:
    mu0 <- par[1]
    #if(mu0 < 0)
    #  mu0 <- 1e-8
    Q <- as.matrix(par[2])
    if(Q[1,1] < 0)
      Q <- as.matrix(0)
    b <- as.matrix(par[3])
    #if(b > 0)
    #  b <- as.matrix(0)
    theta <- par[4]
    # End reading parameters
    cat("theta=",theta,"mu0=",mu0, "Q=",Q, "b=",b,"\n")
    
    res <<- .Call("complik_discrete", dat, dims[1], dims[2], mu0, Q, b, theta, k)
    res <- -1*res
    res
  }
  #maxlik(dat, parameters)
  ## Optimization:
  #optim(par = parameters, 
  #      fn=maxlik, dat = as.matrix(dat), control = list(fnscale=-1, trace=T, ndeps=c(1e-9,1e-9,1e-9,1e-2), maxit=10000), 
  #      method="L-BFGS-B")
  result <- nlm(p = parameters, 
      f=maxlik, dat = as.matrix(dat), gradtol = 1e-26, 
      stepmax = 1, iterlim=10000, steptol = 1e-26, 
      ndigit = 16)
  result
}

