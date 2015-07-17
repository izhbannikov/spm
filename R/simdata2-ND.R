# 2D simulation script
simdata2_ND <- function(N=5, 
                        u=c(3.6461530856,0.2005311687, 10.3993562073, 0.3658859232, 1.1636556317), 
                        R= c(0.9550290288423, -0.0026978184479,  2.365200699e-03, -2.968857693e-05, -0.009103776168, 
                             -0.0129440857139,  0.9754119735264, -2.039885799e-03,  1.077985893e-01,  0.025906794362,
                             0.0146482999044, -0.0074792942366,  9.685838436e-01, -7.320690926e-02, -0.042607607658,
                             0.0003115070409, -0.0008795683512,  1.969713264e-04,  9.916396461e-01, -0.002502486202,
                             0.0050019130410,  0.0057847616518,  9.285817785e-06,  4.981558476e-03,  0.968220352325), 
                        b=c(3.5197233477,9.7709878147, 10.3457493062, 0.5688062855, 4.6505840867), 
                        mumu=1.747957160e-03, 
                        bbb=c(-2.560654665e-05,  1.379145597e-06, -4.796292392e-06,  1.558897881e-06, -3.284264297e-06), 
                        Q=c(1.997835933e-07,  7.500174437e-10,  1.234081917e-08, -4.685963111e-07,  3.773961121e-08,
                            7.500174437e-10,  1.713366327e-09, -8.974449194e-10,-4.560323761e-08, -1.641051801e-09,
                            1.234081917e-08,  -8.974449194e-10,  7.606328447e-09,  1.958545097e-08, -4.984033651e-09,
                            -4.685963111e-07, -4.560323761e-08,  1.958545097e-08, 5.678193990e-07,  3.064454995e-08,
                            3.773961121e-08, -1.641051801e-09, -4.984033651e-09, 3.064454995e-08, 6.965948397e-09), 
                        theta=7.920000000e-02, 
                        tstart=45, 
                        ystart=c(85.63412671,81.74419604,223.4612069, 26.0134818, 50.57496147), 
                        dt=1, 
                        nj=85, 
                        tmax=105,
                        data_names=c("dbp1", "dbp2", "gluc1", "gluc2", "chol1", "chol2"),
                        k=5) {
  
  
  u<-matrix(u,nrow=k,ncol=1,byrow=T)
  b<-matrix(b,nrow=k,ncol=1,byrow=T)
  bbb<-matrix(bbb,nrow=k,ncol=1,byrow=T)
  ystart<-matrix(ystart,nrow=k,ncol=1,byrow=T)
  R <- matrix(R,nrow=k,ncol=k,byrow=T)
  Q <- matrix(Q,nrow=k,ncol=k,byrow=T)
  
  data <- matrix(nrow=0, ncol=(4+length(data_names)))
  colnames(data) <- c("id", "xi", "t1", "t2", data_names)
  
  mu <- function(t,y) {
    mu <- ( mumu + t(y) %*% bbb +t(y) %*% Q %*% y )*exp(theta*t)
    mu
  }
  
  id <- 0
  for(i in 1:N) {
    #i=1
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
          eps[ii,1] <- rnorm(n=1,mean=0,sd=b[ii,1])
        }
        
        y2 <- u + R %*% y1 + eps
      
        new_person <- F
      } else {
        xi <- 1
        y2 <- matrix(nrow=k,ncol=1,NA)
        new_person <- T
      }
      #print(y2)
      dd <- c(y1[1,1], y2[1,1])
      if(k>1) {
        for(ii in 2:k) {
          dd <- c(dd, y1[ii,1], y2[ii,1])
        }
      }
      
      data <- rbind(data, c(id, xi, t1, t2, dd))
      #data <- rbind(data, c(id, xi, t1, t2, y1[1,1], y2[1,1], y1[2,1], y2[2,1]))
      
      if (new_person == F) {
        t1 <- t2
        t2 <- t1 + dt
        
        if(t2 > tmax+1) break;
        
        y1 <- y2
        
        
        if(j == nj) {
          new_person <- T
        }
        j <- j+1
      }
    }
  }
  
  #data <- data[2:dim(data)[1],]
  data
}


# Mortality:
#mumu <- mu0 + f^2*Q
#bbb <- -2*f*Q

# Dynamics:
#u <- -a*f1
#R <- 1+a

#simdata2_ND(N=50,k=3)
