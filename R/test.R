# Test spm
source("~/Projects/spm/R/spmND.R")

#source("Z:/iz12/Projects/spm/rscripts/simdata2-ND.R")

sasdat <- read.sas7bdat("~/Projects/spm/data/fhsm1mortality.sas7bdat", debug=FALSE)
#sasdat <- read.sas7bdat("~/Dropbox/spm/fhsm1mortality.sas7bdat", debug=FALSE)
dat <- sasdat[,1] #pid
dat <- cbind(dat, sasdat[,5]) #sta
dat <- cbind(dat, sasdat[,3]) #tt1
dat <- cbind(dat, sasdat[,4]) #tt3
# 1st dimension:
dat <- cbind(dat, sasdat[,7]) #dbp
dat[2:dim(dat)[1],5] <- dat[1:(dim(dat)[1]-1),5]
dat <- cbind(dat, sasdat[,7]) #dbp
# 2nd dimension:
dat <- cbind(dat, sasdat[,10]) #gluc
dat[2:dim(dat)[1],7] <- dat[1:(dim(dat)[1]-1),7]
dat <- cbind(dat, sasdat[,10]) #gluc
# 3rd dimension:
dat <- cbind(dat, sasdat[,9]) #chol
dat[2:dim(dat)[1],9] <- dat[1:(dim(dat)[1]-1),9]
dat <- cbind(dat, sasdat[,9]) #chol
## 4th dimension:
dat <- cbind(dat, sasdat[,8]) #bmi
dat[2:dim(dat)[1],11] <- dat[1:(dim(dat)[1]-1),11]
dat <- cbind(dat, sasdat[,8]) #bmi
# 5th dimension:
dat <- cbind(dat, sasdat[,6]) #pp
dat[2:dim(dat)[1],13] <- dat[1:(dim(dat)[1]-1),13]
dat <- cbind(dat, sasdat[,6]) #pp


averages = matrix(nrow=1,ncol=5)
averages[1,1] = dat[1,5]
averages[1,2] = dat[1,7]
averages[1,3] = dat[1,9]
averages[1,4] = dat[1,11]
averages[1,5] = dat[1,13]

pid=dat[1,1]
starttime = c(dat[1,3])
for(i in 1:dim(dat)[1]) {
  if(dat[i,1] != pid) {
    dat[i,5] = dat[i,6]
    dat[i,7] = dat[i,8]
    dat[i,9] = dat[i,10]
    dat[i,11] = dat[i,12]
    dat[i,13] = dat[i,14]
    pid = dat[i,1]
    averages<-rbind(averages,c(dat[i,5],dat[i,7], dat[i,9], dat[i,11], dat[i,13]))
    starttime <- c(starttime, dat[i,3])
  }
  if(dat[i,2] > 1) {
    dat[i,2] <- 1
  }
}

avg_dbp=mean(averages[,1])
avg_dbp
avg_gluc=mean(averages[,2],na.rm=T)
avg_gluc
mean(averages[,3],na.rm=T)
mean(averages[,4],na.rm=T)
mean(averages[,5],na.rm=T)
tstart <- mean(starttime)

#dat <- simdata2_ND(N=500,ystart=c(avg_dbp,avg_gluc,230, ),tstart=tstart, data_names = c())
spm(dat, 5)
