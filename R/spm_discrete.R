#'Discrete multi-dimensional optimization
#'@param dat A data table.
#'@param k A number of dimensions.
#'@param theta_range A range of theta parameter (axe displacement of Gompertz function), default: from 0.001 to 0.09 with step of 0.001.
#'@param tol A tolerance threshold for matrix inversion (NULL by default).
#'@return A list of two elements: (1) estimated parameters u, R, b, Sigma, Q, mu0, theta and
#'(2) estimated parameters a, f1, Q, f, b, mu0, theta. Note: b and mu0 from first list are different 
#'from b and mu0 from the second list.
#'@details This function is way much faster that continuous \code{spm_continuous_MD(...)} (but less precise) and used mainly in 
#'estimation a starting point for the \code{spm_continuous_MD(...)}.
#'@examples
#'library(spm)
#'# Reading longitudinal data
#'longdat <- read.csv(system.file("data","longdat.csv",package="spm"))
#'# Prepare data for optimization
#'vitstat <- read.csv(system.file("data","vitstat.csv",package="spm"))
#'data <- prepare_data(longdat=longdat, vitstat=vitstat,interval=1, col.status="IsDead", col.id="ID", col.age="Age", col.age.event="LSmort", covariates=c("DBP"), verbose=T)
#'# Parameters estimation
#'pars <- spm_discrete(data[[2]], k=1, theta_range=seq(0.001,0.09,by=0.001), tol=NULL)
#'pars
spm_discrete <- function(dat,k=1, theta_range=seq(0.001,0.09,by=0.001), tol=NULL, verbose=FALSE) {
  options(digits=10)
  # Logistic regression:
  total_cols <- (1 + k + (k*(k+1))/2) + 2
  result <- matrix(nrow=0, ncol=total_cols,0)
  for(theta in theta_range) {
    ethetat <- exp(theta*dat[,3])
    newdat <- dat[,2] # Outcome
    newdat <- cbind(newdat,ethetat) #x0
    cnames <- c("xi", "x0")
    #cnames <- c("xi")
    
    # A loop for the bt coefficients:
    index_i <- 1
    for(i in seq(1,(k*2-1),2)) {
      newdat <- cbind(newdat,dat[,(4+i)]*ethetat) #x1
      cnames <- c(cnames,paste("x1","_",index_i,sep=''))
      index_i <- index_i + 1
    }
    
    index_i <- 1
    index_j <- 1
    for(i in seq(1,(k*2-1),2)) {
      for(j in seq(i,(k*2-1),2)) {
        newdat <- cbind(newdat,dat[,(4+i)]*dat[,(4+j)]*ethetat)
        cnames <- c(cnames,paste("x2","_", index_i,"_",index_j,sep=''))
        index_j <- index_j + 1
      }
      index_i <- index_i + 1
      index_j <- index_i
    }
    
    colnames(newdat) <- cnames
    
    #reg_formula <- paste("1 -", cnames[1],"~", paste(cnames[2:length(cnames)],collapse='+'))
    reg_formula <- paste("1 -", cnames[1],"~", paste(cnames[2:length(cnames)],collapse='+'), "- 1")
    
    
    res.pois <- glm(reg_formula, data=as.data.frame(newdat), family = poisson(link=log), 
                    control=list(maxit = 250, trace=verbose))
    res <- glm(reg_formula, data=as.data.frame(newdat), family = binomial(link = log), 
               start=coef(res.pois), 
               control=list(maxit = 250, trace=verbose))
    
    coef <- -1*res$coefficients 
    result <- rbind(result, c(theta, coef, logLik(res)[1]))
    
  }

  parameters_glm <- result[which(result[,total_cols] == max(result[,total_cols])),]
  names(parameters_glm) <- c("theta", cnames[2:length(cnames)], "Log Lik")

  # Least-square:
  index_i <- 1
  index_j <- 1
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
  
  #Output parameters:
  parameters_glm <- unname(parameters_glm);
  theta <- parameters_glm[1]
  mu0 <- parameters_glm[2]
  b <- parameters_glm[3:(2+k)]
  #------Q-matrix------#
  Q <- matrix(nrow=k, ncol=k, 0)
  names(Q) <- NULL
  start <- 2+k+1
  end <- start + k*(k+1)/2
  utarray <- parameters_glm[start:end]
  start <- 1
  end <- k
  i <- 1
  for(ii in 1:k) {
    Q[ii,(i:k)] <- utarray[start:end]
    start <- start + k - i + 1
    end <- start + k - i - 1
    i <- i + 1
  }
  for(i in 1:k) {
    for(j in i:k) {
      Q[j,i] <- Q[i,j]
    }
  }
  #------end of calculating of Q-matrix--------#
  parameters_lsq <- unname(parameters_lsq);
  u <- parameters_lsq[,1]
  R <- parameters_lsq[,2:(2+k-1)]
  Sigma <- parameters_lsq[,(2+k)]
  
  pars1 <- list(theta=theta, mu0=mu0, b=b, Q=Q, u=u, R=R, Sigma=Sigma)
  
  # Making a new parameter set:
  QH <- Q
  aH <- R - diag(k)
  bH <- as.matrix(Sigma, nrow=1)
  f1H <- (-1)*u %*% solve(aH, tol=tol)
  fH <- -0.5 * b %*% solve(QH, tol=tol)
  mu0H <- mu0 - fH %*% QH %*% t(fH)
  thetaH <- theta
  pars2 <- list(a=aH, f1=f1H, Q=QH, f=fH, b=bH, mu0=mu0H, theta=thetaH)
  
  pars <- list(pars1=pars1, pars2=pars2)
  invisible(pars)
}

