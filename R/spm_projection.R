spm_projection <- function(x, N, ystart=80, tstart=30, tend=105, dt=1, f=NULL) {
  
  res <- list()
  
  if(length(x) == 1) {
    # Data simulation for time-dependent model
    
    formulas.work <- list(at="-0.05", f1t="80", Qt="2e-8", ft="80", bt="5", mu0t="2e-5")
    
    if (!is.null(f)) {
      for(item in f) {
        formulas.work[[item]] <- formulas[[item]]
      }
    }
    
    #Simulate (project) data:
    res.time_dep <- simdata_time_dep(N=N,f=f,
                                    step=dt/20, tstart=tstart, tend=tend, ystart=ystart, sd0=4, k=1)
    #Compute summary statistics:
    stat <- data.frame(mean=mean(res.time_dep[,5], na.rm=TRUE), 
                        sd=mean(res.time_dep[,5], na.rm=TRUE))
    
    res <- list(time_dependent=res.time_dep, stat=stat)
    
  } else if(length(x) == 2) {
    
    if(dim(x[["Ya2007"]]$a)[1] != length(ystart)) {
      stop("Number of dimensions does not match with the number of values provided in ystart.")
    }
    
    x <- p.discr.model
    N=2000
    tstart=30
    tend=165
    ystart=c(80,25)
    dt=1
    
    # Data simulation for discrete and continuous models
    res.cont <- simdata_cont(N=N, a=x[["Ya2007"]]$a, f1=x[["Ya2007"]]$f1, Q=x[["Ya2007"]]$Q, 
                  f=x[["Ya2007"]]$f, b=x[["Ya2007"]]$b, mu0=x[["Ya2007"]]$mu0, theta=x[["Ya2007"]]$theta,
                  step=dt/20, tstart=tstart, tend=tend, ystart=ystart, sd0=4, k=length(ystart))
    
    res.discr <- simdata_discr(N=N, a=x[["Ya2007"]]$a, f1=x[["Ya2007"]]$f1, Q=x[["Ya2007"]]$Q, 
                 f=x[["Ya2007"]]$f, b=x[["Ya2007"]]$b, mu0=x[["Ya2007"]]$mu0, theta=x[["Ya2007"]]$theta, 
                 ystart=ystart, tstart=tstart, tend=tend, dt=1, k=length(ystart))
    
    stat <- list()
    stat[["continuous"]] <- list()
    for(i in seq(5,length(colnames(res.cont)),by=2)) {
      print(colnames(res.cont)[i])
      stat[["continuous"]][[colnames(res.cont)[i]]][["mean"]] <- mean(res.cont[,i], na.rm = TRUE)
      stat[["continuous"]][[colnames(res.cont)[i]]][["sd"]] <- sd(res.cont[,i], na.rm = TRUE)
    }
    stat[["discrete"]] <- list()
    for(i in seq(5,length(colnames(res.discr)),by=2)) {
      print(colnames(res.discr)[i])
      stat[["discrete"]][[colnames(res.discr)[i]]][["mean"]] <- mean(res.discr[,i], na.rm = TRUE)
      stat[["discrete"]][[colnames(res.discr)[i]]][["sd"]] <- sd(res.discr[,i], na.rm = TRUE)
    }
    
    res <- list(discrete=res.discr, continuous=res.cont, stat=stat)
  }
  
  invisible(res)
}

plot_surv <- function(x) {
  srv <- c()
  n <- max(x[,1])
  for(i in round(min(x[,3]),0):round(max(x[,3]),0)) {
    srv <- c(srv, 
              length(which(x[which(round(x[,3],0) == i), 2] == 0)) / n)
  
  }
  plot(round(min(x[,3]),0):round(max(x[,3]),0), srv, type = "l")
  
}
  