library(stats)

spm_integral <- function(dat,parameters) {
  iteration <- 0
  
  LL <- function(dat, par) {
    a <- par[1] 
    f1 <- par[2]
    Q <- par[3]
    b <- par[4]
    f <- par[5]
    mu0 <- par[6]
    theta <- par[7]
    
    dims <- dim(dat)
    res <- .Call("complik", dat, dims[1], dims[2], a, f1, Q, b, f, mu0, theta)
    
    print(paste("L=",res))
    print(paste("Iter:", iteration, 
                "a=",par[1],
                "f1=",par[2],
                "Q=",par[3],
                "b=",par[4],
                "f=",par[5],
                "mu0=",par[6], 
                "theta=", par[7]))
    
    iteration <<- iteration + 1
  
    res
  }

  # Optimization:
  result <- optim(par = parameters, 
                fn=LL, dat = as.matrix(dat), control = list(fnscale=-1, trace=T, maxit=10000), 
                method="L-BFGS-B", lower = c(-0.1,0,2e-15,0.5,0,2e-15, 0)) 
                #upper = c(0,100,1e-5,10,100,0.5, 1))

  result
}


