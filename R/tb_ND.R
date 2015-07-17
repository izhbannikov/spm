# Testbench

#source("~/Dropbox/spm/simdata2-ND.R")
source("Z:/iz12/Projects/spm/rscripts/simdata2-ND.R")

options(digits=10)

N <- 10000
numsim <- 500
k=5
ests <- matrix(nrow=numsim,ncol=63,0) # Parameter esitmations
colnames(ests) <- c("theta", "mu0", "y1_b1", "y2_b2", "y3_b3", "y4_b4", "y5_b5",
                    "y1_y1", "y1_y2", "y1_y3", "y1_y4", "y1_y5", "y2_y2", "y2_y3", 
                    "y2_y4", "y2_y5", "y3_y3", "y3_y4", "y3_y5", "y4_y4", "y4_y5",
                    "y5_y5", "Log Lik Theta", 
                    "u1", "R11", "R12", "R13", "R14", "R15", "b1", "Log Lik 1", 
                    "u2", "R21", "R22", "R23", "R24", "R25", "b2", "Log Lik 2",
                    "u3", "R31", "R32", "R33", "R34", "R35", "b3", "Log Lik 3",
                    "u4", "R41", "R42", "R43", "R44", "R45", "b4", "Log Lik 4",
                    "u5", "R51", "R52", "R53", "R54", "R55", "b5", "Log Lik 5")

for(step in 1:numsim) {
  # Data simulation:
  dat <- simdata2_ND(N=N,k=5, data_names = c("dbp1","dbp2","gluc1","gluc2","chol1","chol2","bmi1","bmi2","pp1","pp2"), ystart=c(85.63412671,81.74419604,223.4612069, 26.0134818, 50.57496147))
  
  # Logistic regression:
  total_cols <- (1 + k + (k*(k+1))/2) + 2
  #total_cols <- (1 + k + k*k) + 2
  result <- matrix(nrow=0, ncol=total_cols)
  for(theta in seq(0.05,0.1,by=0.001)) {
    newdat <- dat[,2] # Outcome
    newdat <- cbind(newdat,exp(theta*dat[,3])) #x0 - time t1
    
    cnames <- c("xi", "u0")
    # Cycle for bt coefficients:
    index_i <- 1
    for(i in seq(1,(k*2-1),2)) {
      newdat <- cbind(newdat,dat[,(4+i)]*exp(theta*dat[,3])) #y1_b1t
      cnames <- c(cnames,paste("y",index_i,"_", "b",index_i,sep=''))
      index_i <- index_i + 1
    }
    
    index_i <- 1
    index_j <- 1
    for(i in seq(1,(k*2-1),2)) {
      for(j in seq(i,(k*2-1),2)) {
        newdat <- cbind(newdat,dat[,(4+i)]*dat[,(4+j)]*exp(theta*dat[,3]))
        cnames <- c(cnames,paste("y",index_i,"_", "y",index_j,sep=''))
        index_j <- index_j + 1
      }
      index_i <- index_i + 1
      index_j <- index_i
    }
    
    reg_formula <- paste(cnames[1],"~", paste(cnames[2:length(cnames)],collapse='+'), "- 1")
    
    colnames(newdat) <- cnames
    
    res <- glm(reg_formula , data=as.data.frame(newdat)) # -1 means no intercept
    #print(theta)
    #print(summary(res))
    coef <- res$coefficients 
    result <- rbind(result, c(theta, coef, logLik(res)[1]))
  }
  
  parameters_glm <- result[which(result[,total_cols] == max(result[,total_cols])),]
  names(parameters_glm) <- c("theta", cnames[2:length(cnames)], "Log Lik")
  parameters_glm
  
  
  
  # Least-square:
  index_i <- 1
  parameters_lsq <- matrix(nrow=0,ncol=(k+3))
  for(i in seq(1,(k*2-1),2)) {
    newdat2 <- dat[,(5+i)]
    cnames <- c(paste("y2","_",index_i,sep=''))
    for(j in seq(1,(k*2-1),2)) {
      newdat2 <- cbind(newdat2,dat[,(4+j)]) 
      cnames <- c(cnames,paste("y1","_", index_j,sep=''))
      index_j <- index_j  + 1
    }
    colnames(newdat2) <- cnames
    
    reg_formula <- paste(cnames[1],"~", paste(cnames[2:length(cnames)],collapse='+'))
    
    res <- lm(as.formula(reg_formula),data=as.data.frame(newdat2))
    coef <- res$coefficients # intercept, y1
    parameters_lsq <- rbind(parameters_lsq, c(coef, sd(res$residuals), logLik(res)[1]))
    
    
    index_j <- 1
    index_i <- index_i + 1
    
  }
  
  parameters_lsq
  
  
  ests[step, ] <- c(parameters_glm, parameters_lsq[1,], parameters_lsq[2,], parameters_lsq[3,], parameters_lsq[4,], parameters_lsq[5,])
}

head(ests)


true_par <- c(7.920000000e-02,  1.747957160e-03, -2.560654665e-05,  1.379145597e-06, -4.796292392e-06,  1.558897881e-06, -3.284264297e-06,
              1.997835933e-07,  7.500174437e-10,  1.234081917e-08, -4.685963111e-07,  3.773961121e-08,  1.713366327e-09, -8.974449194e-10,
              -4.560323761e-08, -1.641051801e-09,  7.606328447e-09,  1.958545097e-08, -4.984033651e-09,  5.678193990e-07,  3.064454995e-08,
              6.965948397e-09,
              3.6461530856, 0.9550290288423, -0.0026978184479,  2.365200699e-03, -2.968857693e-05, -0.009103776168, 3.5197233477,
              0.2005311687, -0.0129440857139,  0.9754119735264, -2.039885799e-03,  1.077985893e-01,  0.025906794362, 9.7709878147,
              10.3993562073, 0.0146482999044, -0.0074792942366,  9.685838436e-01, -7.320690926e-02, -0.042607607658, 10.3457493062,
              0.3658859232, 0.0003115070409, -0.0008795683512,  1.969713264e-04,  9.916396461e-01, -0.002502486202, 0.5688062855,
              1.1636556317, 0.0050019130410,  0.0057847616518,  9.285817785e-06,  4.981558476e-03,  0.968220352325, 4.6505840867) # theta, mu0, bbb, Q, u0, R, b

#png("~/Dropbox/spm/tb_results_5D_N500_numsim100.png",height=1080, width=1920)

png("/iz12/Projects/spm/rscripts/tb_results_5D_N2000_numsim100.png",height=2*1080, width=2*1920)


par(mfrow=c(6,10))

jj <- 1
for(i in 1:dim(ests)[2]) {
  if(length(grep('Log',colnames(ests)[i])) == 0) {
    ans <- hist(ests[,i],plot=F)
    ymax <- max(ans$density) + 1/10*max(ans$density)
    hist(ests[,i], freq=FALSE, col="green", xlab=colnames(ests)[i], main = paste("Distribution of",colnames(ests)[i]),ylim=c(0, ymax))
    curve(dnorm(x, mean=mean(ests[,i]), sd=sd(ests[,i])), add=TRUE, col="darkblue", lwd=2) 
    abline(v = true_par[jj], col="red", lwd=3)
    text(true_par[jj], y = ymax, labels = paste(paste("True",colnames(ests)[i]),":",true_par[jj]),pos=4, cex=1.5)
    jj <- jj+1
  }
}

dev.off()

