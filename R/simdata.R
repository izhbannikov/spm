simdata <-
function(N=10, aH=-0.5, f1H=80, fH=80, bH=10, QH=1e-05, mu0H=2e-05, thetaH=0.08,m0=80, nj=10, step=0.05, trange=c(65,80), yrange = c(60,100), sd0=4) {
  trapezoid <- function(fun, a, b, n=100,t) {
    h <- (b-a)/n
    x <- seq(a,b,by=h)
    y <- fun(idx=x,t=t)
    s <- h*(y[1]/2 + sum(y[2:n]) + y[n+1]/2)
    s
  }
  
  Q <- function(t) {
    Q <- QH*exp(thetaH*t)
    Q
  }
  
  mu0 <- function(t) {
    mu0 <- mu0H*exp(thetaH*t)
    mu0
  }
  
  myfunc <- function(t, y, parms) {
    dy1 <- aH*(y[1] - f1H) - 2*y[2]*Q(t)*(y[1] - fH)
    dy2 <- 2*aH*y[2] + bH^2 - 2*y[2]^2*Q(t)
    list(c(dy1, dy2))
  }
  
  mu <- function(t,idx) {
    mu <- mu0(t) + Q(t)*(m[idx] - fH)^2 + Q(t)*gamma[idx]
    mu
  }
  
  new_person <- T
  data <- matrix(nrow=1,ncol=5,0)
  record <- 1
  for(i in 1:N) {
    #print(i)
    if(new_person == T) {
      t1 <- runif(1,trange[1], trange[2])
      #print(t1)
      t2 <- t1 + 2*runif(1,0,1) + step
      
      #print(t2)
      while(1) {
        y1 <- rnorm(1,mean=m0, sd=sd0)
        if (y1 >= yrange[1] | y1 <= yrange[2]) {
          break
        } 
      }
      
      out <- euler(times = c(t1,t2), y = c(y1, 0), func = myfunc, parms = NULL)
      m <- out[,2]
      gamma <- out[,3]
      #print(m)
      #print(gamma)
      S <- exp(-1*trapezoid(mu, 1,2,n=1,t2))
      #print(S)
      
      xi <- 0
      if (S < runif(1,0,1)) {
        xi <- 0
      } else {
        xi <- 1
      }
      
      if(xi == 0) {
        y2 <- rnorm(1,mean=m[2], sd=sqrt(gamma))
        new_person <- F
      } else {
        y2 <- NA
        new_person <- T
      }
      
      data <- rbind(data, c(xi, t1, t2, y1, y2))
      record <- record + 1
    } 
    
    j <- 1
    while(new_person == F){
      t1 <- t2
      t2 <- t1 + 2*runif(1,0,1) + step
      y1 <- y2
      
      out <- euler(times = c(t1,t2), y = c(y1, 0), func = myfunc, parms = NULL)
      m <- out[,2]
      gamma <- out[,3]
      S <- exp(-1*trapezoid(mu, 1,2,n=1,t2))
      
      xi <- 0
      if (S < runif(1,0,1)) {
        xi <- 0
      } else {
        xi <- 1
      }
      
      if(xi == 0) {
        y2 <- rnorm(1,mean=m[2], sd=sqrt(gamma))
        new_person <- F
      } else {
        y2 <- NA
        new_person <- T
      }
      
      if(runif(1,0,1) < 0.95) {
        data <- rbind(data, c(xi, data[record,3], t2, data[record,5], y2))
        record <- record + 1
      }
      
      j <- j + 1
      if(j == nj) {
        new_person <- T
      }
    }
    
  }
  
  # One last step:
  data <- data[2:dim(data)[1],]
  data
}
