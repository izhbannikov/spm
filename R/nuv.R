nuv <- function(t, nut,vt, b=1, bb=2e-7, u=4, R=0.95, theta=0.08) {
  
  sigma <- bb #var(rnorm(N, mean=0, sd=bb*exp(theta*t)))
  
  vstar <- 1/(1/vt + bb*exp(theta*t))
  nustar <- nut - vstar*(b*exp(theta*t) + bb*exp(theta*t)*nut)
  
  nu <- u + R*nustar
  
  v <- sigma + R*vstar*R
  
  print(paste(t,u*exp(theta*t),R,sigma,nu,nustar,v,vstar,bb, bb*exp(theta*t)))
  
  c(nu, nustar, v, vstar)
}

set.seed(125)

tstart <- 30
tmax <- 105
nu0 <- 80
v0 <- 0
nustar <- nu0
vstar <- 2e-7
nuvdat <- matrix(nrow=tmax, ncol=5)
colnames(nuvdat) <- c("year","nu","nustar", "v", "vstar")

nut <- nu0
vt <- v0
nuvdat[1,] <- c(tstart,nut, nustar,vt,vstar)
t <- tstart
for(i in 2:tmax) {
  t <- t + 1
  nuvt <- nuv(t, nut,vt)
  nut<-nuvt[1]
  vt <- nuvt[3]
  nuvdat[i,] <- c(t,nut, nuvt[2], vt, nuvt[4])
}

nuvdat

data <- simdata2(1000)
# Error calculation:
diff <- matrix(nrow=1, ncol=5)
for(i in tstart:tmax) {
  avg <- mean(data[which(data[,2]==0 & data[,3]==i),5])
  if (length(which(data[,2]==0 & data[,3]==i) > 1)) {
    
    se <- sd(data[which(data[,2]==0 & data[,3]==i),5])/sqrt(length(which(data[,2]==0 & data[,3]==i)))
  
  } else if(length(which(data[,2]==0 & data[,3]==i) == 1)) {
      
    se <- 0
  
  } else {
    print(paste(i, data[which(data[,2]==0 & data[,3]==i),5]))
    
    break
  }
  
  print(c(i,avg,se))
  diff <- rbind(diff,c(i,avg,se, nuvdat[i,3], (avg-nuvdat[i,3])/se))
}
diff <- diff[2:dim(diff)[1],]
colnames(diff) <- c("age", "mean(y1)", "se(y1)", "nu_star", "(mean(y)-nu_star)/se(y)")
diff
