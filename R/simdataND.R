# ND simulation function
# It uses a, f1, Q, f, b, mu0, theta
sim <- function(N=100, a, f1, Q, f, b, mu0, theta, ystart, tstart=30, tend=105, dt=1, k=1) {
  u_ <- -1*a %*% f1
  R_ <- 1 + a
  epsilon_ <- b
  mu0_ <- mu0 + t(f) %*% Q %*% f
  b_ <- -1*t(f) %*% Q - Q %*% f
  Q_ <- Q
  theta_ <- theta
  
  simulated <- sumdata(N=N, u=u_, R=R_, epsilon=epsilon_, mu0=mu0_, b=b_, Q=Q_, theta=theta_, tstart=tstart, ystart=ystart, dt=dt, tmax=tend, k=k)
  simulated
}

# Function that simulates data using u, R, epsilon, mu0, b, Q, theta
simdata <- function(N=100, # Number of individuals
                        u, 
                        R, 
                        epsilon, 
                        mu0, 
                        b, 
                        Q, 
                        theta, 
                        tstart=30, 
                        ystart, 
                        dt=1, 
                        tmax=105,
                        k=1) {
  
  
  u<-matrix(u,nrow=k,ncol=1,byrow=T)
  epsilon<-matrix(epsilon,nrow=k,ncol=1,byrow=T)
  b<-matrix(b,nrow=k,ncol=1,byrow=T)
  ystart<-matrix(ystart,nrow=k,ncol=1,byrow=T)
  
  data_names <- c()
  for(n in 1:k) {
    data_names <- c(data_names, paste("par",n, "_1", sep=''), paste("par",n, "_2", sep=''))
  }
  
  data <- matrix(nrow=0, ncol=(4+length(data_names)))
  colnames(data) <- c("id", "xi", "t1", "t2", data_names)
  
  #q <- function(t,Q) {
  #  res <- Q %*% exp(theta*t)
  #  res
  #}
  
  mu <- function(t,y) {
    mu <- ( mu0 + t(y) %*% b +t(y) %*% Q %*% y )*exp(theta*t)
    mu
  }
  
  id <- 0
  for(i in 1:N) {
    t1 <- tstart
    t2 <- t1 + dt
    y1 <- ystart
    new_person <- F
    id <- id + 1
    j <- 1
    while(new_person == F) {
      S <- exp(-1*mu(t1,y1))
      
      if(S > runif(1,0,1)) {
        xi <- 0
        eps <- matrix(nrow=k,ncol=1,0)
        for(ii in 1:k) {
          eps[ii,1] <- rnorm(n=1,mean=0,sd=epsilon[ii,1])
        }
        
        y2 <- u + R %*% y1 + eps
      
        new_person <- F
      } else {
        xi <- 1
        y2 <- matrix(nrow=k,ncol=1,NA)
        new_person <- T
      }
      
      dd <- c(y1[1,1], y2[1,1])
      if(k>1) {
        for(ii in 2:k) {
          dd <- c(dd, y1[ii,1], y2[ii,1])
        }
      }
      
      data <- rbind(data, c(id, xi, t1, t2, dd))
      
      if (new_person == F) {
        t1 <- t2
        t2 <- t1 + dt
        
        if(t2 > tmax+1) {
          new_person <- T
          break;
        }
        
        y1 <- y2
        
        #if(j == nj) {
        #  new_person <- T
        #}
        j <- j+1
      }
    }
  }
  
  data
}