#### Multiple imputation with SPM ####

sigma_sq <- function(t1, t2, b) {
  # t2 = t_{j}, t1 = t_{j-1}
  ans <- b %*% (t2-t1)
  ans
}

m <- function(y, t1, t2, a, f1) {
  # y = y_{j-1}, t1 = t_{j-1}, t2 = t_{j}
  ans <- y + a %*% t(y - f1) %*% (t2 - t1)
  ans
}

mu <- function(y,t, mu0, Q, f) {
  ans <- mu0 + (y - f) %*% Q %*% t(y - f)
  ans
}


getNextY <- function(y1, t1, t2, b, a, f1) {
  y2 <- rnorm(length(y1),mean=m(y1, t1, t2, a, f1), sd=sqrt(sigma_sq(t1, t2, b)))
  y2
}

getNextY.discr <- function(y1, u, R, b) {
  eps <- rnorm(1, mean=0.0, sd=b)
  y2 <- u + R %*% y1 + eps;
  y2
}

#'Multiple Data Imputation with SPM
#'@param x A path to the longitudinal dataset with missing observations
#'@param minp Number of imputations. Default: 5
#'@param theta_range A range of parameter theta used for optimization, default: seq(0.01, 0.15, by=0.001).
#'@return A list(imputed, imputations, spm.par)
#'@return imputed An imputed dataset.
#'@return imputations Temporary imputed datasets used in multiple imputaitons.
#'@return spm.par A set of SPM parameters used to in imputation of each of 
#'temporary imputed dataset.
#'@examples \dontrun{
#'library(stpm) 
#'##Data preparation ##
#'data <- simdata_discr(N=1000, dt = 2)
#'miss.id <- sample(x=dim(data)[1], size=round(dim(data)[1]/4)) # ~25% missing data
#'incomplete.data <- data
#'incomplete.data[miss.id,5] <- NA
#'incomplete.data[miss.id,6] <- NA
#'## End of data preparation ##
#'
#'# Estimate parameters from the complete dataset #
#'p <- spm_discrete(data, theta_range = seq(0.075, 0.09, by=0.001))
#'p
#'
#'##### Multiple imputation with SPM #####
#'imp.data <- spm.impute(x=incomplete.data, minp=5, theta_range=seq(0.075, 0.09, by=0.001))$imputed
#'
#'## Estimate SPM parameters from imputed data and compare them to the p ##
#'pp.test <- spm_discrete(imp.data, theta_range = seq(0.075, 0.09, by=0.001))
#'pp.test
#'
#'}
spm.impute <- function(x, minp=5, theta_range=seq(0.01, 0.15, by=0.001)) {
  
  ##### Rough parameter estimations #####
  pp <- spm_discrete(x, theta_range = theta_range)

  ######### Data imputation begins here ############
  datasets <- list() # Keeps imputed datasets
  #datasets.tmp <- list() # Not used
  parameters <- list() # Keeps parameters for each imputation

  # Data copy:
  incomplete.data.1 <- x

  ### First pre-processing stage: impute everything which is possible with current data ###
  for(i in 2:(dim(incomplete.data.1)[1]-1)) {
    ### ID
    id.cur <- incomplete.data.1[i,1]
    id.prev <- incomplete.data.1[i-1,1]
    id.next <- incomplete.data.1[i+1,1]
    ### Row
    row.cur <- incomplete.data.1[i, ]
    row.prev <- incomplete.data.1[i-1, ]
    row.next <- incomplete.data.1[i+1, ]
  
    if( (id.cur == id.prev) & (id.cur == id.next)) { # Somewhere in the middle of person
      for(j in seq(5,length(row.cur), by=2)) {
        if(is.na(row.cur[j])) {
          if(!is.na(row.prev[j+1])) { row.cur[j] <- row.prev[j+1] }
        }
      }
    
      for(j in seq(6,length(row.cur),by=2)) {
        if(is.na(row.cur[j])) {
          if(!is.na(row.next[j-1])) { row.cur[j] <- row.next[j-1] }
        }
      }
    }
    
    if( (id.cur == id.prev) & (id.cur != id.next) ) {
      for(j in seq(5,length(row.cur),by=2)) {
        if(is.na(row.cur[j])) {
          if(!is.na(row.prev[j+1])) { row.cur[j] <- row.prev[j+1] }
        }
      }
    }
    
    if( (id.cur != id.prev) & (id.cur == id.next) ) {
      for(j in seq(6,length(row.cur),by=2)) {
        if(is.na(row.cur[j])) {
          if(!is.na(row.next[j-1])) {
            row.cur[j] <- row.next[j-1]
          }
        }
      }
    }
  
    
    incomplete.data.1[i, ] <- row.cur
  }
  

  ### Second pre-processing stage: impute the rest (prediction by simulation) ###
  for(inp in 1:minp) {
    incomplete.data.1.tmp <- incomplete.data.1
    
    #### First, let us handle the first row ####
    row.cur <- incomplete.data.1.tmp[1, ]
    row.next <- incomplete.data.1.tmp[2, ]
    t1 <- row.cur[3]; t2 <- row.cur[4]
    
    for(j in seq(5,length(row.cur),by=2)) {
      if(is.na(row.cur[j])) {
        meanY <- mean(incomplete.data.1[, j], na.rm = TRUE)
        y.start <- rnorm(1, mean = meanY, pp$Ya2007$b[((j-5) %/% 2)+1])
        row.cur[j] <- y.start
      }
    
      if(any(is.na(row.cur[seq(6, length(row.cur),by=2)]))) {
        y1 <- row.cur[seq(5,length(row.cur),by=2)]
        y.next <- getNextY(y1, t1, t2, pp$Ya2007$b, pp$Ya2007$a, pp$Ya2007$f1)
        row.cur[which(is.na(row.cur))] <- y.next[(which(is.na(row.cur)) - 6 ) %/% 2 + 1]
      }
      for(j in seq(6, length(row.cur),by=2)) {
        if(is.na(row.next[j-1])) {row.next[j-1] <- row.cur[j]}
      }
    }
    
    incomplete.data.1.tmp[1, ] <- row.cur
    incomplete.data.1.tmp[2, ] <- row.next

    for(i in 2:(dim(incomplete.data.1.tmp)[1]-1)) {
      ### ID
      id.cur <- incomplete.data.1.tmp[i,1]
      id.prev <- incomplete.data.1.tmp[i-1,1]
      id.next <- incomplete.data.1.tmp[i+1,1]
      ### Row
      row.cur <- incomplete.data.1.tmp[i, ]
      row.prev <- incomplete.data.1.tmp[i-1, ]
      row.next <- incomplete.data.1.tmp[i+1, ]
      t1 <- row.cur[3]; t2 <- row.cur[4]
      
      #####
      if( (id.cur == id.prev) & (id.cur == id.next)) { # Somewhere in the middle of person
        if(any(is.na(row.cur[seq(6, length(row.cur),by=2)]))) {
          y1 <- row.cur[seq(5,length(row.cur),by=2)]
          y.next <- getNextY(y1, t1, t2, pp$Ya2007$b, pp$Ya2007$a, pp$Ya2007$f1)
          row.cur[which(is.na(row.cur))] <- y.next[(which(is.na(row.cur)) - 6 ) %/% 2 + 1]
        }
        for(j in seq(6, length(row.cur),by=2)) {
          if(is.na(row.next[j-1])) {row.next[j-1] <- row.cur[j]}
        }
      }
      #####
      if( (id.cur != id.prev) ) {
        for(j in seq(5,length(row.cur),by=2)) {
          if(is.na(row.cur[j])) {
            meanY <- mean(incomplete.data.1[, j], na.rm = TRUE)
            y.start <- rnorm(1, mean = meanY, pp$Ya2007$b[((j-5) %/% 2)+1])
            row.cur[j] <- y.start
          }
        }
    
        if(any(is.na(row.cur[seq(6, length(row.cur),by=2)])) & row.cur[2] == 0) {
          y1 <- row.cur[seq(5,length(row.cur),by=2)]
          y.next <- getNextY(y1, t1, t2, pp$Ya2007$b, pp$Ya2007$a, pp$Ya2007$f1)
          row.cur[which(is.na(row.cur))] <- y.next[(which(is.na(row.cur)) - 6 ) %/% 2 + 1]
        }
        
        for(j in seq(6,length(row.cur),by=2)) {
          if(is.na(row.next[j-1]) & (id.cur == id.next)) { row.next[j-1] <- row.cur[j] }
        }
      }
    
      
      incomplete.data.1.tmp[i, ] <- row.cur
      incomplete.data.1.tmp[i+1, ] <- row.next
    }

    ### Last record in a dataset ###
    Nrec <- dim(incomplete.data.1.tmp)[1]
    id.cur <- incomplete.data.1.tmp[Nrec,1]
    id.prev <- incomplete.data.1.tmp[Nrec-1,1]
    row.cur <- incomplete.data.1.tmp[Nrec, ]
    row.prev <- incomplete.data.1.tmp[Nrec-1, ]
    t1 <- row.cur[3]; t2 <- row.cur[4]
    if(any(is.na(row.cur[seq(6, length(row.cur),by=2)])) & row.cur[2] == 0) {
      y1 <- row.cur[seq(5,length(row.cur),by=2)]
      y.next <- getNextY(y1, t1, t2, pp$Ya2007$b, pp$Ya2007$a, pp$Ya2007$f1)
      row.cur[which(is.na(row.cur))] <- y.next[(which(is.na(row.cur)) - 6 ) %/% 2 + 1]
    }
    
    incomplete.data.1.tmp[Nrec, ] <- row.cur
  
    ### Saving current imputed copy:
    #datasets.tmp[[inp]] <- incomplete.data.1.tmp
    ### Estimate parameters for this datasets:
    pp.1 <- spm_discrete(incomplete.data.1.tmp, theta_range = theta_range)
    ### Saving current set of parameters:
    parameters[[inp]] <- pp.1
    #print(paste(pp.1$Ya2007$a, pp.1$Ya2007$Q, p$Ya2007$a, p$Ya2007$Q))
    
  }

  ##### Filing the missing data base on determined parameters #####
  for(imp in 1:minp) {
    incomplete.data.1.tmp <- incomplete.data.1
    ### Second pre-processing stage: impute the rest (prediction by simulation) ###
    #### First, let us handle the first row ####
    row.cur <- incomplete.data.1.tmp[1, ]
    row.next <- incomplete.data.1.tmp[2, ]
    for(j in seq(5,length(row.cur),by=2)) {
      if(is.na(row.cur[j])) {
        meanY <- mean(incomplete.data.1[, j], na.rm = TRUE)
        y.start <- rnorm(1, mean = meanY, sd=parameters[[minp]]$Ya2007$b[((j-5) %/% 2)+1])
        row.cur[j] <- y.start
      }
    
      if(any(is.na(row.cur[seq(6, length(row.cur),by=2)])) & row.cur[2]==0) {
        y1 <- row.cur[seq(5,length(row.cur),by=2)]
        y.next <- getNextY(y1, t1, t2, pp$Ya2007$b, pp$Ya2007$a, pp$Ya2007$f1)
        row.cur[which(is.na(row.cur))] <- y.next[(which(is.na(row.cur)) - 6 ) %/% 2 + 1]
      }
      for(j in seq(6, length(row.cur),by=2)) {
        if(is.na(row.next[j-1])) {row.next[j-1] <- row.cur[j]}
      }
    }
    
    incomplete.data.1.tmp[1, ] <- row.cur
    incomplete.data.1.tmp[2, ] <- row.next
  
    for(i in 2:(dim(incomplete.data.1.tmp)[1]-1)) {
      ### ID
      id.cur <- incomplete.data.1.tmp[i,1]
      id.prev <- incomplete.data.1.tmp[i-1,1]
      id.next <- incomplete.data.1.tmp[i+1,1]
      ### Row
      row.cur <- incomplete.data.1.tmp[i, ]
      row.prev <- incomplete.data.1.tmp[i-1, ]
      row.next <- incomplete.data.1.tmp[i+1, ]
    
      #####
      if( (id.cur == id.prev) & (id.cur == id.next)) { # Somewhere in the middle of person
        if(any(is.na(row.cur[seq(6, length(row.cur),by=2)]))) {
          t1 <- row.cur[3]; t2 <- row.cur[4]
          y1 <- row.cur[seq(5,length(row.cur),by=2)]
          y.next <- m(y1, t1, t2, pp$Ya2007$a, pp$Ya2007$f1)
          row.cur[which(is.na(row.cur))] <- y.next[(which(is.na(row.cur)) - 6 ) %/% 2 + 1]
        }
        for(j in seq(6, length(row.cur),by=2)) {
          if(is.na(row.next[j-1])) {row.next[j-1] <- row.cur[j]}
        }
        
      }
      #####
      if( (id.cur != id.prev) ) {
        for(j in seq(5,length(row.cur),by=2)) {
          if(is.na(row.cur[j])) {
            meanY <- mean(incomplete.data.1[, j], na.rm = TRUE)
            y.start <- rnorm(1, mean = meanY, pp$Ya2007$b[((j-5) %/% 2)+1])
            row.cur[j] <- y.start
          }
        }
        
        if(any(is.na(row.cur[seq(6, length(row.cur),by=2)])) & row.cur[2] == 0) {
          t1 <- row.cur[3]; t2 <- row.cur[4]
          y1 <- row.cur[seq(5,length(row.cur),by=2)]
          y.next <- m(y1, t1, t2, pp$Ya2007$a, pp$Ya2007$f1)
          row.cur[which(is.na(row.cur))] <- y.next[(which(is.na(row.cur)) - 6 ) %/% 2 + 1]
        }
        
        for(j in seq(6,length(row.cur),by=2)) {
          if(is.na(row.next[j-1]) & (id.cur == id.next)) { row.next[j-1] <- row.cur[j] }
        }
        
        
      }
    
      incomplete.data.1.tmp[i, ] <- row.cur
      incomplete.data.1.tmp[i+1, ] <- row.next
    }
  
    ### Last record in a dataset ###
    Nrec <- dim(incomplete.data.1.tmp)[1]
    id.cur <- incomplete.data.1.tmp[Nrec,1]
    id.prev <- incomplete.data.1.tmp[Nrec-1,1]
    row.cur <- incomplete.data.1.tmp[Nrec, ]
    row.prev <- incomplete.data.1.tmp[Nrec-1, ]
    if(any(is.na(row.cur[seq(6, length(row.cur),by=2)])) & row.cur[2] == 0) {
      y1 <- row.cur[seq(5,length(row.cur),by=2)]
      y.next <- m(y1, t1, t2, pp$Ya2007$a, pp$Ya2007$f1)
      row.cur[which(is.na(row.cur))] <- y.next[(which(is.na(row.cur)) - 6 ) %/% 2 + 1]
    }
    
    incomplete.data.1.tmp[Nrec, ] <- row.cur
    
    ### Saving current imputed copy:
    datasets[[inp]] <- incomplete.data.1.tmp
  }

  ### Summarizing imputed datasets together by averaging missing values, i.e. 'completion' ###

  final.dataset <- incomplete.data
  for(j in seq(5,dim(final.dataset)[2])) {
    data.tmp <- matrix(nrow = dim(final.dataset)[1], ncol=0)
    for(i in 1:minp) {
      data.tmp <- cbind(data.tmp, datasets[[i]][,j])
    }
    data.tmp.2 <- apply(X = data.tmp, FUN = median, MARGIN = 1, na.rm=T)
    final.dataset[,j] <- data.tmp.2
  }
 
  res <- list(imputed=final.dataset, imputations=datasets, spm.par=parameters)
  res
}


