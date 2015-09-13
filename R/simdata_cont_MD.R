#' Simulation function for continuous trait.
#' @param N Number of individuals.
#' @param a A k by k matrix, which characterize the rate of the adaptive response.
#' @param f1 A particular state, which if a deviation from the normal (or optimal). This is a vector with length of k.
#' @param Q A matrix k by k, which is a non-negative-definite symmetric matrix.
#' @param f A vector-function (with length k) of the normal (or optimal) state.
#' @param b A diffusion coefficient, k by k matrix.
#' @param mu0 mortality at start period of time.
#' @param theta A displacement coefficient of the Gompertz function.
#' @param ystart A vector with length equal to number of dimensions used, defines starting values of covariates.
#' @param tstart A number that defines starting time (30 by default).
#' @param tend A number, defines final time (105 by default).
#' @param k number of dimensions (k = 1 by default).
#' @return A table with simulated data.
#' @examples
#' library(spm)
#' dat <- simdata_cont(N=2500)
#' dat
#
simdata_cont_MD <- function(N=10, aH=-0.05, f1H=80, QH=2e-07, fH=80, bH=5, mu0H=2e-05, thetaH=0.08,
                         step=0.05, tstart=30, tend=105, ystart=80, sd0=4, k=1) {
    
  aH<-matrix(aH,nrow=k,ncol=k,byrow=T)
  f1H<-matrix(f1H,nrow=k,ncol=1,byrow=T)
  QH<-matrix(QH,nrow=k,ncol=k,byrow=T)
  fH<-matrix(fH,nrow=k,ncol=1,byrow=T)
  bH<-matrix(bH,nrow=k,ncol=1,byrow=T)
  ystart<-matrix(ystart,nrow=k,ncol=1,byrow=T)
  
  Q <- function(t) {
    Q <- QH*exp(thetaH*t)
    Q
  }
    
  mu <- function(t, par) {
    #par[1] - m, par[2] - gamma
    hfH <- fH - par[[1]]
    hf1H <- f1H - par[[1]]
      
    mu0Ht <- mu0H*exp(thetaH*t)
    QH_gamma1 <- QH %*% par[[2]]
    mu <- mu0Ht + t(hfH) %*% QH %*% hfH + sum(diag(QH_gamma1))
    mu
  }
  
  
  func1 <- function(t, y) {
    # y[[1]] = m, y[[2]] = gamma
    hfH <- fH - y[[1]]
    hf1H <- f1H - y[[1]]
    dm <- -1.0*aH%*%hf1H + 2.0*y[[2]]%*%Q(t)%*%hfH
    dgamma <- aH%*%y[[2]] + y[[2]]%*%t(aH) + bH%*%t(bH) - 2.0*y[[2]]%*%Q(t)%*%y[[2]] 
    list(dm, dgamma)
  }
  
  new_person <- T
  data <- matrix(nrow=1,ncol=(4+2*k),NA)
  record <- 1
  id <- 1
  for(i in 1:N) {
    if(new_person == T) {
      t1 <- runif(1,tstart, tend) # Starting time
      t2 <- t1 + 2*runif(1,0,1) + step # Ending time
        
      # Starting point
      y1 <- matrix(unlist(lapply(1:k, function(n) {rnorm(1,mean=ystart[n,1], sd=sd0)} )), 
                   nrow=k, ncol=1, byrow=T)
      nsteps <- 2.00
      h=(t2-t1)/nsteps
        
      # Integration:
      s <- h/3*(-1)*mu(t1, list(y1,matrix(0,nrow=k,ncol=k)))
        
      t <- t1
      out <- list(y1,matrix(0,nrow=k,ncol=k))
      yfin <- list()
      ytmp <- list()
      
      for(j in 1:nsteps) {
        # Runge-Kutta method:
        k1ar <- func1(t,out)
        yfin[[1]] <- out[[1]] + h/6.00*k1ar[[1]]
        yfin[[2]] <- out[[2]] + h/6.00*k1ar[[2]]
        ytmp[[1]] <- out[[1]] + h/2.00*k1ar[[1]]
        ytmp[[2]] <- out[[2]] + h/2.00*k1ar[[2]]
        
        k2ar <- func1(t,ytmp)
        yfin[[1]] <- yfin[[1]] + h/3.00*k2ar[[1]]
        yfin[[2]] <- yfin[[2]] + h/3.00*k2ar[[2]]
        ytmp[[1]] <- out[[1]] + h/2.00*k2ar[[1]]
        ytmp[[2]] <- out[[2]] + h/2.00*k2ar[[2]]
        
        k3ar <- func1(t,ytmp)
        yfin[[1]] <- yfin[[1]] + h/3.00*k3ar[[1]]
        yfin[[2]] <- yfin[[2]] + h/3.00*k3ar[[2]]
        ytmp[[1]] <- out[[1]] + h*k3ar[[1]]
        ytmp[[2]] <- out[[2]] + h*k3ar[[2]]
        
        k4ar <- func1(t,ytmp)
        out[[1]] <- yfin[[1]] + h/6.00*k4ar[[1]]
        out[[2]] <- yfin[[2]] + h/6.00*k4ar[[2]]
        
        t <- t + h
          
        # Integration:
        if (j == nsteps) {
          ifactor <- 1
        } else {
          if ((j %% 2) == 0) {
            ifactor <- 2
          } else {
            ifactor <- 4
          }
        }
        s <- s + ifactor*h/3.00*(-1)*mu(t,out)
          
      }
        
      m <- out[[1]]
      gamma <- out[[2]]
        
      S <- exp(s)
        
      xi <- 0
      if (S > runif(1,0,1)) {
        xi <- 0
      } else {
        xi <- 1
      }
      
      if(xi == 0) {
        y2 <- matrix(unlist(lapply(1:k, function(n) {rnorm(1,mean=m[n,1], sd=sqrt(gamma[n,n]))} )), 
               nrow=k, ncol=1, byrow=T)
        new_person <- F
      } else {
        y2 <- matrix(nrow=k, ncol=1, NA)
        new_person <- T
      }
      
      cov <- unlist(lapply(1:k, function(n) {c(y1[n,1], y2[n,1])}))
      data <- rbind(data, c(id, xi, t1, t2, cov))
      record <- record + 1
        
      if(new_person == T) {
          id <- id + 1
        }
      } 
      
      while(new_person == F){
        t1 <- t2
        t2 <- t1 + 2*runif(1,0,1) + step
        y1 <- y2
        
        nsteps <- 2
        h=(t2-t1)/nsteps
        
        # Integration:
        s <- h/3*(-1)*mu(t1, list(y1,matrix(0,nrow=k,ncol=k)))
        
        t <- t1
        out <- list(y1,matrix(0,nrow=k,ncol=k))
        for(j in 1:nsteps) {
          # Runge-Kutta method:
          k1ar <- func1(t,out)
          yfin[[1]] <- out[[1]] + h/6.00*k1ar[[1]]
          yfin[[2]] <- out[[2]] + h/6.00*k1ar[[2]]
          ytmp[[1]] <- out[[1]] + h/2.00*k1ar[[1]]
          ytmp[[2]] <- out[[2]] + h/2.00*k1ar[[2]]
          
          k2ar <- func1(t,ytmp)
          yfin[[1]] <- yfin[[1]] + h/3.00*k2ar[[1]]
          yfin[[2]] <- yfin[[2]] + h/3.00*k2ar[[2]]
          ytmp[[1]] <- out[[1]] + h/2.00*k2ar[[1]]
          ytmp[[2]] <- out[[2]] + h/2.00*k2ar[[2]]
          
          k3ar <- func1(t,ytmp)
          yfin[[1]] <- yfin[[1]] + h/3.00*k3ar[[1]]
          yfin[[2]] <- yfin[[2]] + h/3.00*k3ar[[2]]
          ytmp[[1]] <- out[[1]] + h*k3ar[[1]]
          ytmp[[2]] <- out[[2]] + h*k3ar[[2]]
          
          k4ar <- func1(t,ytmp)
          out[[1]] <- yfin[[1]] + h/6.00*k4ar[[1]]
          out[[2]] <- yfin[[2]] + h/6.00*k4ar[[2]]
          
          t <- t + h
          
          # Integration:
          if (j == nsteps) {
            ifactor <- 1
          } else {
            if ((j %% 2) == 0) {
              ifactor <- 2
            } else {
              ifactor <- 4
            }
          }
          s <- s + ifactor*h/3.00*(-1)*mu(t,out)
          
        }
        
        m <- out[[1]]
        gamma <- out[[2]]
        
        S <- exp(s)
        
        xi <- 0
        if (S > runif(1,0,1)) {
          xi <- 0
          y2 <- matrix(unlist(lapply(1:k, function(n) {rnorm(1,mean=m[n,1], sd=sqrt(gamma[n,n]))} )), 
                       nrow=k, ncol=1, byrow=T)
          new_person <- F
          cov <- unlist(lapply(seq(1,k), function(n) {c(y1[n,1], y2[n,1])}))
          data <- rbind(data, c(id, xi, data[record,4], t2, cov))
          record <- record + 1
          
        } else {
          xi <- 1
          y2 <- matrix(nrow=k, ncol=1, NA)
          new_person <- T
          cov <- unlist(lapply(1:k, function(n) {c(y1[n,1], y2[n,1])}))
          data <- rbind(data, c(id, xi, data[record,4], t2, cov))
          record <- record + 1
          id <- id + 1
        }
        
        if(t2 > tend & new_person == F) {
          new_person <- T
          id <- id + 1
        }
        
      }
      
    }
    
    # One last step:
    data <- data[2:dim(data)[1],]
    colnames(data) <- c("id","xi","t1","t2", unlist(lapply(1:k, function(n) {c(paste("y_", n, "_1", sep=""), paste("y_", n, "_2",sep="") )} )) )
    invisible(data)
}