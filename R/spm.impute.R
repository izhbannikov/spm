#' An internal function to compute sigma square analytically
#' @param t1 t1
#' @param t2 t2
#' @param b b (see Yashin et. al, 2007)
#' @return sigma_square (see Akushevich et. al, 2005)
sigma_sq <- function(t1, t2, b) {
  # t2 = t_{j}, t1 = t_{j-1}
  ans <- b %*% (t2-t1)
  ans
}

#' An internal function to compute m from 
#' @param y Current value of Y
#' @param t1 t1
#' @param t2 t2
#' @param a a (see Yashin et. al, 2007)
#' @param f1 f1 (see Yashin et. al, 2007)
#' @return m m (see Yashin et. al, 2007)
m <- function(y, t1, t2, a, f1) {
  # y = y_{j-1}, t1 = t_{j-1}, t2 = t_{j}
  ans <- y + a %*% t(y - f1) %*% (t2 - t1)
  ans
}

#' An internal function to compute mu
#' @param y Current value of y
#' @param mu0 mu0 (see Yashin et. al, 2007)
#' @param b b (see Yashin et. al, 2007)
#' @param Q Q (see Yashin et. al, 2007)
#' @param theta theta (see Yashin et. al, 2007)
#' @param tt t (time)
#' @return mu Next value of mu
mu <- function(y, mu0, b, Q, theta, tt) {
  ans <- (mu0 + t(y) %*% b + t(y) %*% Q %*% y)*exp(theta*tt)
  ans
}

#mu <- function(tt, y1, gamma1, f, f1, mu0, theta, Q) {
#  hf <- f - y1
#  hf1 <- f1 - y1
#  #if(gomp) {
#  #  mu0Ht = mu0H*exp(thetaH*t);
#  #} else {
#  #  mu0Ht = mu0H;
#  #}
#  #print(t)
#  #print(theta)
#  mu0Ht <- mu0*exp(theta*tt)
#  QH_gamma1 <- Q %*% gamma1
#  mu <- mu0Ht + (t(hf) %*% Q) %*% hf + sum(diag(QH_gamma1))
#  mu
#}

#' An internal function to compute next value of physiological variable Y
#' @param y1 y1
#' @param t1 t1
#' @param t2 t2
#' @param b b (see Yashin et. al, 2007)
#' @param a a (see Yashin et. al, 2007)
#' @param f1 f1 (see Yashin et. al, 2007)
#' @return y.next Next value of y
getNextY.cont2 <- function(y1, t1, t2, b, a, f1) {
  y2 <- rnorm(length(y1),mean=m(y1, t1, t2, a, f1), sd=sqrt(sigma_sq(t1, t2, b)))
  #y2 <- rnorm(length(y1),mean=y1, sd=sqrt(sigma_sq(t1, t2, b)))
  y2
}

#' An internal function to compute the next value of physiological variable Y
#' based on discrete-time model (Akushevich et. al., 2005)
#' @param y1 y1
#' @param u u (see Akushevich et. al, 2005)
#' @param R R (see Akushevich et. al, 2005)
#' @param Sigma Sigma (see Akushevich et. al, 2005)
#' @return y.next Next value of y
getNextY.discr <- function(y1, u, R, Sigma) {
  #eps<-matrix(nrow=dim(R)[1], ncol=1)
  #eps[,1] <- sapply(1:length(Sigma), function(i) {rnorm(1, mean=0.0, sd=Sigma[i])})
  #for(i in 1:length(Sigma)) {
  #  eps[i,1] <- rnorm(1, mean=0.0, sd=Sigma[i])
  #}
  #y2 <- u + R %*% y1  + eps
  y2 <- rnorm(length(y1), mean=getNextY.discr.m(y1, u, R), sd=Sigma)
  y2
}

#' An internal function to compute next m based on dicrete-time model 
#' @param y1 y1
#' @param u u
#' @param R R
#' @return m Next value of m (see Yashin et. al, 2007)
getNextY.discr.m <- function(y1, u, R) {
  m <- u + R %*% y1
  m
}

#' An internal function to compute previous m based on discrete-time model
#' @param y2 y2
#' @param u u
#' @param R R
#' @return m Next value of m (see Yashin et. al, 2007)
getPrevY.discr.m <- function(y2, u, R) {
  m <- solve(R) %*% (y2 - u)
  m
}

#' An internal function to compute previous value of
#' physiological variable Y based on 
#' discrete-time model
#' @param y2 y2
#' @param u u
#' @param R R
#' @param Sigma Sigma
#' @return y1 Previous value of y
getPrevY.discr <- function(y2, u, R, Sigma) {
  eps<-matrix(nrow=dim(R)[1], ncol=1)
  eps[,1] <- sapply(1:dim(eps)[1], function(i) {rnorm(1, mean=0.0, sd=Sigma[i])})
  
  y1 <- solve(R) %*% (y2 - u - eps)
  y1
}


#' An internal function to compute m and gamma based on 
#' continuous-time model (Yashin et. al., 2007)
#' @param tt tt - time
#' @param y y 
#' @param a a (see Yashin et. al, 2007)
#' @param f1 f1 (see Yashin et. al, 2007)
#' @param Q Q (see Yashin et. al, 2007)
#' @param f f (see Yashin et. al, 2007)
#' @param b b (see Yashin et. al, 2007)
#' @param theta theta
#' @return list(m, gamma) Next values of m and gamma (see Yashin et. al, 2007)
func1 <- function(tt, y, a, f1, Q, f, b, theta) {
    #a <- pars[1]; f1 <- pars[2]; Q <- pars[3]; f <- pars[4]; b <- pars[5]; theta <- pars[6];
    hf <- f - y[[1]]
    hf1 <- f1 - y[[1]]
    res <- c()
    dm <- -1.00 * (a %*% hf1) + (2.00 * y[[2]]) %*% Q %*% hf
    dgamma <- a %*% y[[2]] + y[[2]] %*% t(a) + b %*% t(b) - 2.00 * ((y[[2]] %*% Q) %*% y[[2]])
    
    res <- list(m=dm, gamma=dgamma)
    res
}

#' An internal function to compute next Y based on 
#' continous-time model (Yashin et. al., 2007)
#' @param y1 y1
#' @param t1 t1
#' @param t2 t2
#' @param a a (see Yashin et. al, 2007)
#' @param f1 f1 (see Yashin et. al, 2007)
#' @param Q Q (see Yashin et. al, 2007)
#' @param f f (see Yashin et. al, 2007)
#' @param b b (see Yashin et. al, 2007)
#' @param mu0 mu (see Yashin et. al, 2007)
#' @param theta theta (see Yashin et. al, 2007)
#' @param u (see Akushevich et. al, 2007)
#' @param R (see Akushevich et. al, 2007)
#' @return y.next Next value of Y
getNextY.cont <- function(y1, t1, t2, a, f1, Q, f, b, mu0, theta, u, R) {
  nsteps <- 2
  tdiff <- t2-t1
  h <- tdiff/nsteps
  gamma1 <- matrix(nrow=dim(a)[1], ncol=dim(a)[1], 0) # set gamma1 to zero matrix
  #gamma1 <- matrix(nrow=1, ncol=1, 0) # set gamma1 to zero matrix
  tt <- t1
  
  out <- list()
  out[[1]] <- y1
  out[[2]] <- gamma1
  
  
  for(j in 1:nsteps) {
    #Runge-Kutta method:
    
    yn <- out
    k1 <- func1(tt, yn, a, f1, Q, f, b, theta)
    
    yn[[1]] <- out[[1]] + h/2 * k1[[1]]
    yn[[2]] <- out[[2]] + h/2 * k1[[2]]
    k2 <- func1(tt+h/2, yn, a, f1, Q, f, b, theta)
    
    yn[[1]] <- out[[1]] + h/2 * k2[[1]]
    yn[[2]] <- out[[2]] + h/2 * k2[[2]]
    k3 <- func1(tt+h/2, yn, a, f1, Q, f, b, theta)
    
    yn[[1]] <- out[[1]] + h * k3[[1]]
    yn[[2]] <- out[[2]] + h * k3[[2]]
    k4 <- func1(tt+h, yn, a, f1, Q, f, b, theta)
    
    out[[1]] <- out[[1]] + h/6 * (k1[[1]] + 2*k2[[1]] + 2*k3[[1]] + k4[[1]])
    out[[2]] <- out[[2]] + h/6 * (k1[[2]] + 2*k2[[2]] + 2*k3[[2]] + k4[[2]])
    
    tt <- tt + h
  }
  
  m2 <- out[[1]]
  gamma2 <- out[[2]]
  
  # New y2:
  y2 <- matrix(nrow=dim(m2)[1], ncol=1, 0)
  for(ii in 1:dim(m2)[1]) {
    y2[ii,1] <- rnorm(1, mean(m2[ii,1]), sd=sqrt(gamma2[ii,ii])) 
  }
  y2
}

#'Multiple Data Imputation with SPM
#'@param dataset A longitudinal dataset with missing observations
#'@param minp Number of imputations. Default: 5
#'@param theta_range A range of parameter theta used for optimization, default: seq(0.01, 0.15, by=0.001).
#'@return A list(imputed, imputations)
#'@return imputed An imputed dataset.
#'@return imputations Temporary imputed datasets used in multiple imputaitons.
#'@export
#'@examples \dontrun{
#'library(stpm) 
#'##Data preparation ##
#'data <- simdata_discr(N=1000, dt = 2)
#'miss.id <- sample(x=dim(data)[1], size=round(dim(data)[1]/4)) # ~25% missing data
#'incomplete.data <- data
#'incomplete.data[miss.id,5] <- NA
#'incomplete.data[miss.id-1,6] <- NA
#'## End of data preparation ##
#'
#'# Estimate parameters from the complete dataset #
#'p <- spm_discrete(data, theta_range = seq(0.075, 0.09, by=0.001))
#'p
#'
#'##### Multiple imputation with SPM #####
#'imp.data <- spm.impute(dataset=incomplete.data, 
#'                       minp=5, 
#'                       theta_range=seq(0.075, 0.09, by=0.001))$imputed
#'head(imp.data)
#'## Estimate SPM parameters from imputed data and compare them to the p ##
#'pp.test <- spm_discrete(imp.data, theta_range = seq(0.075, 0.09, by=0.001))
#'pp.test
#'}
spm.impute <- function(dat, 
                       col.id=1, 
                       col.status=2,
                       col.age=3, 
                       col.age.event=3, 
                       covariates=4,
                       minp=5, 
                       theta_range=seq(0.01, 0.2, by=0.001), 
                       format="short") 
{
    
    # Check input parameters for correctness
    if(class(dat) != "data.frame")
    {
        stop("Class of dataset must be a 'data.frame'.")
    }
  
    datasets <- list() # To keep imputed datasets
    
    if(format == "short")
    {
        # Prepare data to be in format id xi t1 t2 y y.next
        dataset <- prepare_data_cont(dat, 
                                col.id.ind=col.id, 
                                col.status.ind=col.status,
                                col.age.ind=col.age, 
                                col.age.event.ind=col.age.event, 
                                col.covar.ind=covariates, 
                                dt=1, 
                                verbose=FALSE,
                                impute=FALSE)
    } else if(format == "long")
    {
        cov.tmp <- c(covariates, covariates)
        cov.tmp[seq(1,length(cov.tmp), 2)] <- covariates
        cov.tmp[seq(2,length(cov.tmp), 2)] <- covariates + 1
        dataset <- dat[, c(col.id, col.status, col.age, col.age.event, cov.tmp)]
    } else 
    {
        stop("Format is incorrectly defined.")
    }
    
    # Estimate parameters from raw data
    pp <- spm_discrete(dataset, theta_range = theta_range)
  
    # Individual IDs
    ids <- unique(dataset[,1])
  
    for(m in 1:minp)
    {
        x <- dataset
        
        Ncol <- dim(x)[2]
        for(k in ids) 
        {
            ########## Forward #########
            df <- x[which(x[,1] == k), ]
      
            if(length(df[,1]) == 1) 
            {
                Nrec <- 1
            } else 
            {
                Nrec <- length(df[,1])
            }
      
            if(Nrec == 1) 
            {
                row.cur <- df
                for(j in seq(5,Ncol,by=2)) 
                {
                    if(is.na(row.cur[j]) & !is.na(row.cur[j+1])) 
                    {
                        kkk <- ((j-5) %/% 2)+1
            
                        #y.start <- rnorm(1, mean = mean(x[,j], na.rm = T), sd=pp$Ak2005$Sigma[((j-5) %/% 2)+1])
                        #row.cur[j] <- y.start
            
                        y2 <- df[1,j+1]
                        y1 <- getPrevY.discr.m(t(as.matrix(y2)), pp$dmodel$u[kkk], pp$dmodel$R[kkk,kkk])
                        df[1,j] <- y1
                        row.cur[j] <- y1
                    } else if(is.na(row.cur[j]) & is.na(row.cur[j+1]))
                    {
                        y.start <- rnorm(1, mean = mean(x[,j], na.rm = T), sd=pp$dmodel$Sigma[((j-5) %/% 2)+1])
                        row.cur[j] <- y.start
                    }
                }
                
                if(any(is.na(row.cur[seq(6, Ncol,by=2)])) & row.cur[2] == 0) 
                {
                    y1 <- row.cur[seq(5,Ncol,by=2)]
                    y.next <- getNextY.discr.m(t(as.matrix(y1)), pp$dmodel$u, pp$dmodel$R)
                    #y.next <- getNextY.discr(t(as.matrix(y1)), pp$Ak2005$u, pp$Ak2005$R, pp$Ak2005$Sigma)
                    row.cur[which(is.na(row.cur))] <- y.next[(which(is.na(row.cur)) - 6 ) %/% 2 + 1]
                }
        
                df <- row.cur
                x[which(x[,1] == k), ] <- df
                
                next
            }
      
            # Check the first row #
            row.cur <- df[1, ]
            row.next <- df[2, ]
  
            for(j in seq(5,Ncol,by=2)) 
            {
                if(is.na(row.cur[j])) 
                {
                    kkk <- ((j-5) %/% 2)+1
                    ##meanY <- mean(x[, j], na.rm = TRUE)
                    #medianY <- median(x[, j], na.rm = TRUE)
                    ##y.start <- rnorm(1, mean = meanY, sd=pp$Ak2005$Sigma[((j-5) %/% 2)+1])
                    #y.start <- medianY #rnorm(1, mean = meanY, sd=pp$Ak2005$Sigma[((j-5) %/% 2)+1])
          
                    for(ii in 1:Nrec) {
                        if(!is.na(df[ii,j]))
                            break
                    }
                    
                    if(ii != Nrec) {
                        for(iii in (ii-1):1) {
                            y2 <- df[iii,j+1]
                            y1 <- getPrevY.discr.m(t(as.matrix(y2)), pp$dmodel$u[kkk], pp$dmodel$R[kkk,kkk])
                            df[iii,j] <- y1
                            if(iii != 1) {
                                df[iii-1,j+1] <- y1
                            }
                        }
                    } else {
                        meanY <- mean(x[, j], na.rm = TRUE)
                        y.start <- rnorm(1, mean = meanY, sd=pp$dmodel$Sigma[((j-5) %/% 2)+1])
                        df[1,j] <- y.start
                    }
                } 
              
                if(is.na(row.cur[j+1])) {
                    row.cur[j+1] <- getNextY.discr.m(t(as.matrix(row.cur[j])), pp$dmodel$u, pp$dmodel$R)
                    df[1,j+1] <- row.cur[j+1]
                }
            }
      
            # Check the rest #
            #df[1, ] <- row.cur
            #df[2, ] <- row.next
      
            #### Preprocessing of the rest of df ####
            for(i in 2:Nrec) 
            {
                row.cur <- df[i,]
                row.prev <- df[i-1,]
                if(i != Nrec & Nrec >2) {row.next <- df[i+1,]}
          
                for(j in seq(5, Ncol,by=2)) 
                {
                    if(is.na(row.cur[j])) 
                    {
                        row.cur[j] <- row.prev[j+1]
                    } else 
                    {
                        row.prev[j+1] <- row.cur[j]
                    }
                }
          
                for(j in seq(6, Ncol,by=2)) 
                {
                    if(is.na(row.cur[j])) 
                    {
                        if(i != Nrec & Nrec >2) {row.cur[j] <- row.next[j-1]}
                    } else 
                    {
                        if(i != Nrec & Nrec >2) {row.next[j-1] <- row.cur[j]}
                    }
                }
          
                df[i-1, ] <- row.prev
                df[i, ] <- row.cur
                if(i != Nrec & Nrec >2) {df[i+1, ] <- row.next}
            }
      
            #if(dim(x[which(x[,1] == k), ])[1] != dim(df)[1]) {
            #  print("!!")
            #  print(x[which(x[,1] == k), ])
            #  print(df)
            #  print("???")
            #}
      
        
            #### Main imputation loop ####
        
            for(i in 2:Nrec) 
            {
                row.cur <- df[i, ]
                if(i != Nrec & Nrec >2) {row.next <- df[i+1, ]}
                #print(row.cur)
                y1 <- row.cur[seq(5,Ncol,by=2)]
                if(any(is.na(row.cur[seq(6, Ncol,by=2)]))) 
                {
                    y.next <- getNextY.discr.m(t(as.matrix(y1)), pp$dmodel$u, pp$dmodel$R)
                    #y.next <- getNextY.discr(t(as.matrix(y1)), pp$Ak2005$u, pp$Ak2005$R, pp$Ak2005$Sigma)
                    #y.next <- getNextY.cont(y1, row.cur[3], row.cur[4], pp$Ya2007$a, pp$Ya2007$f1, pp$Ya2007$Q, pp$Ya2007$f, pp$Ya2007$b, pp$Ya2007$mu0, pp$Ya2007$theta, pp$Ak2005$u, pp$Ak2005$R)
                    #y.next <- getNextY.cont(t(as.matrix(y1)), t1=row.cur[3], t2=row.cur[4], a=pp$Ya2007$a, f1=pp$Ya2007$f1, Q=pp$Ya2007$Q, f=pp$Ya2007$f, b=pp$Ya2007$b, mu0=pp$Ya2007$mu0, theta=pp$Ya2007$theta, u=pp$Ak2007$u, R=pp$Ak2005$R) 
                    #y.next <- getNextY.cont2(as.matrix(y1), row.cur[3], row.cur[4], pp$Ya2007$b, pp$Ya2007$a, pp$Ya2007$f1)
                    for(j in seq(6,Ncol,by=2)) 
                    {
                        if(is.na(row.cur[j])) { row.cur[j] <- y.next[(j - 6 ) %/% 2 + 1] }
                        if(i != Nrec & Nrec >2) {if(is.na(row.next[j-1])) { row.next[j-1] <- row.cur[j] }}
                    }
                }
    
                df[i, ] <- row.cur
                if(i != Nrec & Nrec >2) {df[i+1, ] <- row.next}
            }
      
            ### Last record in a dataset ###
            row.cur <- df[Nrec, ]
            if(any(is.na(row.cur[seq(6, Ncol,by=2)])) & row.cur[2] == 0) 
            {
                y1 <- row.cur[seq(5,Ncol,by=2)]
                y.next <- getNextY.discr.m(t(as.matrix(y1)), pp$dmodel$u, pp$dmodel$R)
                #y.next <- getNextY.discr(t(as.matrix(y1)), pp$Ak2005$u, pp$Ak2005$R, pp$Ak2005$Sigma)
                for(j in seq(6, Ncol, by=2)) 
                {
                    if(is.na(row.cur[j])) 
                    {
                        row.cur[j] <- y.next[(j - 6) %/% 2 + 1]
                    }
                }
            }
      
            df[Nrec, ] <- row.cur
      
            tryCatch({
                x[which(x[,1] == k), ] <- df
            },error=function(e) {
                print(e)
                print("x:")
                print(x[which(x[,1] == k), ])
                print("df:")
                print(df)
            }, warning=function(w){
                print(w)
                print("x:")
                print(x[which(x[,1] == k), ])
                print("df:")
                print(df)
                x[which(x[,1] == k), ] <- df
            })
        }
        #################################################################################
        datasets[[m]] <- x
    }
  
    ### Summarizing imputed datasets together by averaging missing values, i.e. 'completion' ###
    final.dataset <- dataset
    for(j in seq(5,dim(final.dataset)[2])) 
    {
        data.tmp <- matrix(nrow = dim(final.dataset)[1], ncol=0)
        for(i in 1:minp) 
        {
            data.tmp <- cbind(data.tmp, datasets[[i]][,j])
        }
        data.tmp.2 <- apply(X = data.tmp, FUN = mean, MARGIN = 1, na.rm=T)
        final.dataset[,j] <- data.tmp.2
    }
  
    res <- list(imputed=final.dataset, imputations=datasets)
    res
}