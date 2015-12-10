library(stats)

#'Continuous one-dimensional optimization
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
#'@details \code{spm_integral_1D} runs much slower that discrete but more precise and can handle time intervals with different lengths.
#'@examples
#'library(spm)
#'# Reading the data:
#'longdat <- read.csv(system.file("data","longdat.csv",package="spm"))
#'vitstat <- read.csv(system.file("data","vitstat.csv",package="spm"))
#'dd <- prepare_data(longdat=longdat, vitstat=vitstat,interval=1, col.status="IsDead", col.id="ID", col.age="Age", col.age.event="LSmort", covariates=c("DBP"), verbose=T)
#'data <- dd[[1]][,2:6]
#'#Parameters estimation:
#'pars <- spm_continuous_1D(dat=data,a=-0.05, f1=80, Q=2e-8, f=80, b=5, mu0=2e-5, theta=0.08)
#'pars
spm_continuous_1D <- function(dat, a=-0.05, f1=80, Q=2e-8, f=80, b=5, mu0=2e-5, theta=0.08) {
  
  res_prev <- NULL
  pars_prev <<- c(a, f1, Q, f, b, mu0, theta)
  iteration <- 0
  
  lower_bound <<- c(-0.5,0,1e-12,1e-6,1e-6,1e-6, 1e-4)
  upper_bound <<- c(0, Inf, 1e-7, Inf, Inf, 1, 0.1)
  
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
    
    res
  }
  
  
  #maxlik(as.matrix(dat), pars_prev)
  # Optimization:
  result <- optim(par = pars_prev, 
                  fn=maxlik, dat = as.matrix(dat), control = list(fnscale=-1, trace=T, maxit=10000, factr=1e-16, ndeps=c(1e-12,1e-12,1e-16,1e-12,1e-12,1e-12,1e-12)), 
                  method="L-BFGS-B", lower = lower_bound, upper = upper_bound)
  
  result
}

