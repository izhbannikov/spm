library(stats)


spm_integral_1D <- function(dat,parameters) {
  
  res_prev <- NULL
  pars_prev <<- parameters
  iteration <- 0
  
  lower_bound <<- c(-1,0,1e-16,0,0,0, 0)
  upper_bound <<- c(0, Inf, 1e-6, Inf, Inf, Inf, 0.1)
  
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
  
  
  maxlik(as.matrix(dat), pars_prev)
  # Optimization:
  result <- optim(par = pars_prev, 
                  fn=maxlik, dat = as.matrix(dat), control = list(fnscale=-1, trace=T, maxit=10000, factr=1e-16, ndeps=c(1e-12,1e-12,1e-16,1e-12,1e-12,1e-12,1e-12)), 
                  method="L-BFGS-B", lower = lower_bound, upper = upper_bound)
  
  result
}

