#' Simulation function for continuous trait with time-dependant coefficients.
#' @param N Number of individuals.
#' @param formulas : a list of formulas that define age (time) - dependency. Default: list(at="a", f1t="f1", Qt="Q*exp(theta*t)", ft="f", bt="b", mu0t="mu0*exp(theta*t)")
#' @param ystart A starting value of covariates.
#' @param tstart A number that defines starting time (30 by default).
#' @param tend A number, defines final time (105 by default).
#' @return A table with simulated data.
#' @examples
#' library(spm)
#' dat <- simdata_time_dep(N=2500)
#' dat
simdata_time_dep <- function(N=10,f=list(at="-0.05", f1t="80", Qt="2e-8", ft="80", bt="5", mu0t="2e-5"),
                         step=0.05, tstart=30, tend=105, ystart=80, sd0=4, k=1) {
  formulas <- f  
  at <- NULL
  f1t <- NULL
  Qt <- NULL
  ft <- NULL
  bt <- NULL
  mu0t <- NULL
  
  comp_func_params <- function(astring, f1string, qstring, fstring, bstring, mu0string) {
    at <<- eval(bquote(function(t) .(parse(text = astring)[[1]])))
    f1t <<- eval(bquote(function(t) .(parse(text = f1string)[[1]]))) 
    Qt <<- eval(bquote(function(t) .(parse(text = qstring)[[1]])))
    ft <<- eval(bquote(function(t) .(parse(text = fstring)[[1]])))
    bt <<- eval(bquote(function(t) .(parse(text = bstring)[[1]])))
    mu0t <<- eval(bquote(function(t) .(parse(text = mu0string)[[1]])))
  }
  
    
  mu <- function(t, par) {
    hfH <- ft(t)-par[[1]];
    hf1H <- f1t(t)-par[[1]];
    
    mu <- mu0t(t)+hfH^2*Qt(t)+Qt(t)*par[[2]]
    mu
  }
  
  func1 <- function(t, y) {
    hfH <- ft(t)-y[[1]]
    hf1H <- f1t(t)-y[[1]]
    dy1 <- -1.0*at(t)*hf1H + 2.0*y[[2]]*Qt(t)*hfH
    dy2 <- 2.0*at(t)*y[[2]] + bt(t) - 2.0*y[[2]]^2*Qt(t) #dy2 <- 2*aH*y[2] + bH^2 - 2*y[2]^2*Q(t)
    list(dy1, dy2)
  }
  
  comp_func_params(formulas$at, formulas$f1t, formulas$Qt, formulas$ft, formulas$bt, formulas$mu0t)
  
  new_person <- T
  data <- matrix(nrow=1,ncol=(4+2*k),NA)
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
      s <- h/3*(-1)*mu(t1, list(y1,0))
        
      t <- t1
      out <- list(y1,0)
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
        #y2 <- matrix(unlist(lapply(1:k, function(n) {rnorm(1,mean=m[n,1], sd=sqrt(gamma[n,n]))} )), 
        #       nrow=k, ncol=1, byrow=T)
        y2 <- rnorm(1,mean=m, sd=sqrt(gamma))
        new_person <- F
      } else {
        y2 <- NA
        new_person <- T
      }
      
      cov <- c(y1, y2)
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
        s <- h/3*(-1)*mu(t1, list(y1,0))
        
        t <- t1
        out <- list(y1,0)
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
          y2 <- rnorm(1,mean=m, sd=gamma)
          new_person <- F
          cov <- c(y1, y2)
          data <- rbind(data, c(id, xi, data[record,4], t2, cov))
          record <- record + 1
          
        } else {
          xi <- 1
          y2 <- NA
          new_person <- T
          cov <- c(y1, y2)
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