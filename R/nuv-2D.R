nuv <- function(t, nut,vt, 
                b=c(3.523178055e+00,9.784105643e+00), 
                bb=c(1.551743641e-07, -1.202172282e-09, 
                     -1.202172282e-09, 6.632778093e-10),  
                u=c(3.765964237e+00,2.642588552e+00), 
                R= c(9.554959987e-01,-4.118460781e-03, 
                     -1.424732045e-03, 9.804356540e-01),
                theta=8.200000000e-02) {
  
  nut=matrix(nut,nrow=2,ncol=1, byrow=F)
  vt=matrix(vt,nrow=2,ncol=2, byrow=F)
  u=matrix(u,nrow=2,ncol=1, byrow=F)
  b=matrix(b,nrow=2,ncol=1, byrow=F)
  bb=matrix(bb,nrow=2,ncol=2,byrow=F)
  R <- matrix(R,nrow=2,ncol=2,byrow=F)
  
  sigma <- bb #var(rnorm(N, mean=0, sd=bb*exp(theta*t)))
  
  if (all(vt == 0)) {
    vstar <- vt
  } else {
    vstar <- solve((solve(vt) + bb*exp(theta*t)))
  }
  
  nustar <- nut - vstar %*% ( b*exp(theta*t) + t(t(nut) %*% (bb*exp(theta*t))) )
  
  nu <- u + R %*% nustar
  
  v <- sigma + R %*% vstar %*% t(R)
  
  #print(paste(t,u*exp(theta*t),R,sigma,nu,nustar,v,vstar,bb, bb*exp(theta*t)))
  
  list(nu, nustar, v, vstar)
}

#source("Z:/iz12/Projects/spm/rscripts/simdata2-2D.R")
source("~/Dropbox/spm/simdata2-2D.R")

tstart <- 44.54951798
tmax <- 100
nu0 <- c(85.63412671,81.74419604)
v0 <- c(0,0,
        0,0)
nustar <- nu0

vstar <- c(1.551743641e-07, -1.202172282e-09, 
           -1.202172282e-09, 6.632778093e-10)
vstar <- matrix(vstar,nrow=2,ncol=2,byrow=F)
nuvdat <- matrix(nrow=tmax, ncol=13)
colnames(nuvdat) <- c("year","nu1", "nu2","nustar1", "nustar2", "v11", "v12", "v21", "v22", "vstar11", "vstar12", "vstar21", "vstar22")

nut <- nu0
vt <- v0
nuvdat[1,] <- c(tstart,nut[1], nut[2], nustar[1], nustar[2],vt[1], vt[2], vt[3], vt[4], vstar[1,1], vstar[1,2], vstar[2,1], vstar[2,2])
t <- tstart
for(i in 2:tmax) {
  t <- t + 1
  nuvt <- nuv(t, nut,vt)
  nut<-nuvt[[1]]
  nustar <- nuvt[[2]]
  vt <- nuvt[[3]]
  vstar <- nuvt[[4]]
  nuvdat[i,] <- c(t,nut[1,1], nut[2,1], nustar[1,1], nustar[2,1], vt[1,1], vt[1,2], vt[2,1], vt[2,2], vstar[1,1], vstar[1,2], vstar[2,1], vstar[2,2])
}

nuvdat

data <- simdata2_2D(500,tstart=tstart,ystart=c(85.63412671,81.74419604))
# Error calculation:
diff <- matrix(nrow=1, ncol=9)
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
  
  # For the second dimension:
  avg2 <- mean(data[which(data[,2]==0 & data[,3]==i),7])
  if (length(which(data[,2]==0 & data[,3]==i) > 1)) {
    
    se2 <- sd(data[which(data[,2]==0 & data[,3]==i),7])/sqrt(length(which(data[,2]==0 & data[,3]==i)))
    
  } else if(length(which(data[,2]==0 & data[,3]==i) == 1)) {
    
    se2 <- 0
    
  } else {
    print(paste(i, data[which(data[,2]==0 & data[,3]==i),7]))
    
    break
  }
  
  #print(c(i,avg,se))
  diff <- rbind(diff,c(i,avg,se, nuvdat[(i-tstart+1),4], (avg-nuvdat[(i-tstart+1),4])/se, avg2,se2, nuvdat[(i-tstart+1),5], (avg2-nuvdat[(i-tstart+1),5])/se2))
}
diff <- diff[2:dim(diff)[1],]
colnames(diff) <- c("age", "mean(dpb)", "se(dpb)", "nu_star_dbp", "(mean(dbp)-nu_star_dbp)/se(dbp)",  "mean(gluc)", "se(gluc)", "nu_star_gluc", "(mean(gluc)-nu_star_gluc)/se(gluc)")
diff
