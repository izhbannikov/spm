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
simdata_cont <- function(N=10, aH=-0.05, f1H=80, QH=2e-07, fH=80, bH=5, mu0H=2e-05, thetaH=0.08,
                         step=0.05, tstart=30, tend=105, ystart=80, sd0=4) {
     
    Q <- function(t) {
      Q <- QH*exp(thetaH*t)
      Q
    }
    
    mu <- function(t, par) {
      hfH <- fH-par[1];
      hf1H <- f1H-par[1];
      
      mu0Ht <- mu0H*exp(thetaH*t);
      mu <- mu0Ht+hfH^2*QH+QH*par[2];
      mu
    }
    
    func1 <- function(t, y) {
      hfH <- fH-y[1];
      hf1H <- f1H-y[1];
      dy1 <- -1.0*aH*hf1H + 2.0*y[2]*Q(t)*hfH
      dy2 <- 2.0*aH*y[2] + bH - 2.0*y[2]^2*Q(t) #dy2 <- 2*aH*y[2] + bH^2 - 2*y[2]^2*Q(t)
      c(dy1, dy2)
    }
    
    new_person <- T
    data <- matrix(nrow=1,ncol=6,NA)
    record <- 1
    id <- 1
    for(i in 1:N) {
      if(new_person == T) {
        t1 <- runif(1,tstart, tend) # Starting time
        t2 <- t1 + 2*runif(1,0,1) + step # Ending time
        
        # Starting point
        y1 <- rnorm(1,mean=ystart, sd=sd0)
        
        nsteps <- 2.00
        h=(t2-t1)/nsteps
        
        # Integration:
        s <- h/3*(-1)*mu(t1,c(y1,0))
        
        t <- t1
        out <- c(y1,0)
        for(j in 1:nsteps) {
          # Runge-Kutta method:
          k1ar <- func1(t,out)
          yfin <- out + h/6.00*k1ar
          ytmp <- out + h/2.00*k1ar
          
          k2ar <- func1(t,ytmp)
          yfin <- yfin + h/3.00*k2ar
          ytmp <- out + h/2.00*k2ar
          
          k3ar <- func1(t,ytmp)
          yfin <- yfin + h/3.00*k3ar
          ytmp <- out + h*k3ar
          
          k4ar <- func1(t,ytmp)
          out <- yfin + h/6.00*k4ar
          
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
          s <- s + ifactor*h/3.00*(-1)*mu(t,c(out[1],out[2]))
          
        }
        
        m <- out[1]
        gamma <- out[2]
        
        S <- exp(s)
        
        xi <- 0
        if (S > runif(1,0,1)) {
          xi <- 0
        } else {
          xi <- 1
        }
        
        if(xi == 0) {
          y2 <- rnorm(1,mean=m, sd=sqrt(gamma))
          new_person <- F
        } else {
          y2 <- NA
          new_person <- T
        }
        
        data <- rbind(data, c(id, xi, t1, t2, y1, y2))
        record <- record + 1
        
        if(new_person == T) {
          id <- id + 1
        }
      } 
      
      jj <- 1
      while(new_person == F){
        t1 <- t2
        t2 <- t1 + 2*runif(1,0,1) + step
        y1 <- y2
        
        nsteps <- 2.00
        h=(t2-t1)/nsteps
        
        # Integration:
        s <- h/3*(-1)*mu(t1,c(y1,0))
        
        t <- t1
        out <- c(y1,0)
        for(j in 1:nsteps) {
          # Runge-Kutta method:
          k1ar <- func1(t,out)
          yfin <- out + h/6.00*k1ar
          ytmp <- out + h/2.00*k1ar
          
          k2ar <- func1(t,ytmp)
          yfin <- yfin + h/3.00*k2ar
          ytmp <- out + h/2.00*k2ar
          
          k3ar <- func1(t,ytmp)
          yfin <- yfin + h/3.00*k3ar
          ytmp <- out + h*k3ar
          
          k4ar <- func1(t,ytmp)
          out <- yfin + h/6.00*k4ar
          
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
          s <- s + ifactor*h/3.00*(-1)*mu(t,c(out[1],out[2]))
          
        }
        
        m <- out[1]
        gamma <- out[2]
        
        S <- exp(s)
        
        xi <- 0
        if (S > runif(1,0,1)) {
          xi <- 0
          y2 <- rnorm(1,mean=m, sd=sqrt(gamma))
          new_person <- F
        } else {
          xi <- 1
          y2 <- NA
          new_person <- T
          id <- id + 1
        }
        
        if(t2 > tend) {
          new_person <- T
          id <- id + 1
        } else {
          data <- rbind(data, c(id, xi, data[record,4], t2, data[record,6], y2))
          record <- record + 1
        } 
        
      }
      
    }
    
    # One last step:
    data <- data[2:dim(data)[1],]
    colnames(data) <- c("id","xi","t1","t2","y1","y2")
    invisible(data)
}