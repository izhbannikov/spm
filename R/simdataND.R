#' Multi-dimension simulation function
#' It uses a, f1, Q, f, b, mu0 and theta
#' as input parameters.
#' @param N Number of individuals
#' @param k number of dimensions (k = 1 by default).
#' @param a A k by k matrix, which characterize the rate of the adaptive response.
#' @param f1 A particular state, which is a deviation from the normal (or optimal). This is a vector with length of k.
#' @param Q A matrix k by k, which is a non-negative-definite symmetric matrix.
#' @param f A vector-function (with length k) of the normal (or optimal) state.
#' @param b A diffusion coefficient, k by k matrix.
#' @param mu0 mortality at start period of time.
#' @param theta A displacement coefficient of the Gompertz function.
#' @param ystart A vector with length equal to number of dimensions used, defines starting values of covariates.
#' @param tstart A number that defines starting time (30 by default).
#' @param tend A number, defines final time (105 by default).
#' @param dt A time step (1 by default).
#' @return A table with simulated data.
#' @examples
#' library(spm)
#' data <- sim(N=1000, ystart=c(75, 94), k=1)
#' head(data)
sim_discrete <- function(N=100, a=-0.05, f1=80, Q=2e-8, f=80, b=5, mu0=1e-5, theta=0.08, ystart=c(80), tstart=30, tend=105, dt=1, k=1) {
  # Re-calculating parameters:
  u_ <- matrix((f1 %*% (-1*a)), nrow=k, byrow=T)
  #matrix(u,nrow=k,ncol=1,byrow=T)
  R_ <- matrix((diag(k) + a), nrow=k, ncol=k, byrow=T)
  epsilon_ <- matrix(b, nrow=k, ncol=1, byrow=T)
  mu0_ <- mu0 + (f %*% Q) %*% t(f)
  b_ <- matrix((-1*f %*% Q - f %*% Q), nrow=k, ncol=1, byrow=T)
  Q_ <- matrix(Q, nrow=k, ncol=k, byrow=T)
  theta_ <- theta
  ystart = matrix(ystart, nrow=k, ncol=1)
  simulated <- simdata(N=N, u=u_, R=R_, epsilon=epsilon_, mu0=mu0_, b=b_, Q=Q_, theta=theta_, tstart=tstart, ystart=ystart, dt=dt, tmax=tend, k=k)
  #simulated = .Call("simdata_ND", N, u_, R_, epsilon_, mu0_, b_, Q_, theta_, tstart, ystart, tend, k);
  
  data_names <- c()
  for(n in 1:k) {
    data_names <- c(data_names, paste("par",n, "_1", sep=''), paste("par",n, "_2", sep=''))
  }
  colnames(simulated) <- c("id", "xi", "t1", "t2", data_names)
  
  invisible(simulated)
}

#' Function that simulates data using u, R, epsilon, mu0, b, Q, theta
#' @param N Number of individuals
#' @param k Number of dimensions (k = 1 by default).
#' @param u A drift vector with length of k.
#' @param R A k by k regression matrix.
#' @param epsilon A time-dependent normally distributed random vector (size=k).
#' @param mu0 mortality at start period of time.
#' @param b A diffusion coefficient, k by k matrix.
#' @param Q A matrix k by k, which is a non-negative-definite symmetric matrix.
#' @param theta A displacement coefficient of the Gompertz function.
#' @param ystart A vector with length equal to number of dimensions used, defines starting values of covariates.
#' @param tstart A number that defines starting time (30 by default).
#' @param tend A number, defines final time (105 by default).
#' @param dt A time step (1 by default).
#' @return A table with simulated data.
#' @examples
#' library(spm)
#' data <- simdata(N=1000, ystart=c(75, 94), k=1)
#' head(data)
simdata <- function(N=100, # Number of individuals
                        u=8, 
                        R=0.95, 
                        epsilon=5, 
                        mu0=2e-5, 
                        b=10, 
                        Q=2e-8, 
                        theta=0.08, 
                        tstart=30, 
                        ystart=80, 
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