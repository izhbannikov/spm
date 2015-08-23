library(stats)

spm_integral <- function(dat,parameters) {
  iteration <- 0
  LL <- function(dat, par) {
    #print(dat)
    aH <- par[1] 
    f1H <- par[2]
    QH <- par[3]
    bH <- par[4]
    fH <- par[5]
    mu0H <- par[6]
    #thetaH <- par[7]
    thetaH_ <- 0.08
  
    dims <- dim(dat)
    res <- .Call("complik", dat, dims[1], dims[2], aH, f1H, QH, bH, fH, mu0H, thetaH_)
    
    print(paste("L=",res))
    print(paste("Iter:", iteration, "aH=",par[1],"f1H=",par[2],"QH=",par[3],"bH=",par[4],"fH=",par[5],"mu0H=",par[6]))
    iteration <<- iteration + 1
  
    res
  }


  ## Optimization:
  #iteration <<- 0
  #
  result <- optim(par = parameters, 
                fn=LL, dat = as.matrix(dat), control = list(fnscale=-1, trace=T, maxit=10000), 
                method="L-BFGS-B", lower = c(-0.1,0,2e-15,0.5,0,2e-15), 
                upper = c(0,100,1e-5,10,100,0.5))



  
  #ans=LL(dat, parameters)
  result
}


