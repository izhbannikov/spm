# Preparing data for stochastic process model
prepare_data <- function(longdat, vitstat, interval=1, col.status="IsDead", col.id="ID", col.age="Age", col.age.event="LSmort", covariates=c("DBP", "BMI", "DBP1", "DBP2", "Weight", "Height"), verbose=T) {
  
  # Parsing input parameters in order to check for errors:
  if( !(col.status %in% colnames(vitstat)) ) {
    stop(paste("Status column",col.status, "not found in vitstat table. Aborting."))
  }
  if( !(col.id %in% colnames(vitstat) || col.id %in% colnames(longdat)) ) {
    stop(paste("ID column",col.id, "not found in vitstat and/or longdat tables. Aborting."))
  }
  if( !(col.age %in% colnames(longdat)) ) {
    stop(paste("Age column",col.age, "not found in longdat table. Aborting."))
  }
  if( !(col.age.event %in% colnames(vitstat)) ) {
    stop(paste("Event column",col.age.event, "not found in vitstat table. Aborting."))
  }
  for(c in covariates) {
    if( !(c %in% colnames(longdat)) ) {
      stop(paste("Covariate",c, "not found. Aborting."))
    }
  }
  
  #-----------Done parsing imput parameters---------------------#
  
  fill_last <- function(x) {
    na_idx <- which(is.na(x))
    unique_elements <- unique(x[-na_idx])
    set_diff <- unique_elements[length(unique_elements)]
    x[na_idx] <- set_diff
    x
  }
  
  approx2p <- function(t1, y1, t2, y2, t) {
    # 2-point linear interpolation:
    y <- y1 + (t - t1)/(t2 - t1)*(y2 -y1)
    y
  }
  
  # Interpolation
  dt <- interval
  tt <- matrix(nrow=0, ncol=4)
  par <- matrix(nrow=0, ncol=length(covariates))
  
  splitted <- split(longdat, longdat[[col.id]])
  vitstat.splitted <- split(vitstat, vitstat[[col.id]])
  
  for(iii in 1:length(splitted)) {
    if(!is.na(vitstat.splitted[[iii]][[col.age.event]]) & !is.na(vitstat.splitted[[iii]][[col.status]]) ) {
      if(verbose)
        print(iii)
    
      id <- splitted[[iii]][[col.id]][1]
      nrows <- (tail(splitted[[iii]][[col.age]], n=1) - splitted[[iii]][[col.age]][1])/dt + 1
      # Perform approximation:
      t1.approx <- matrix(ncol=4, nrow=nrows)
      t1.approx[,1] <- id
      t1.approx[,2] <- 0
      t1.approx[nrows,2] <- vitstat.splitted[[iii]][[col.status]][1] #Last value
      t1.approx[,3] <- seq(splitted[[iii]][[col.age]][1], splitted[[iii]][[col.age]][length(splitted[[iii]][[col.age]])], by=dt)
      if(nrows > 1) {
        t1.approx[,4] <- c(t1.approx[,3][2:nrows], vitstat.splitted[[iii]][[col.age.event]][1])
      } else {
        t1.approx[,4] <- vitstat.splitted[[iii]][[col.age.event]][1]
      }
    
      tt <- rbind(tt,t1.approx)
      par1.approx <- matrix(ncol=length(covariates), nrow=nrows, NA)
    
      j <- 1
      for(name in covariates) {
        name <- covariates[j]
        if ( (length(splitted[[iii]][[name]]) > 1) & (length(which(!is.na(splitted[[iii]][[name]]))) > 0) ) {
          if(length(which(!is.na(splitted[[iii]][[name]]))) == 1) {
            splitted[[iii]][[name]] <- fill_last(splitted[[iii]][[name]])
          }
          # Fill NAs by linear approximation with approx():
          nn <- length(splitted[[iii]][[name]])
          splitted[[iii]][[name]] <- approx(splitted[[iii]][[name]],n=nn)$y
          par1.approx[,j] <-  approx(splitted[[iii]][[name]], n=nrows)$y
        }
      
        j <- j + 1
      
      }
      par <- rbind(par,par1.approx)
    }
  }
  
  ans=cbind(tt,par)
  colnames(ans) <- c("ID", "CASE", "T1", "T3", covariates)
  #print(head(ans))
  ans <- ans[rowSums( matrix(is.na(ans[,5:dim(ans)[2]]), ncol=length(covariates),byrow=T)) !=length(covariates),]
  
  ans_final <- ans
  if(length(which(is.na(ans[,5:dim(ans)[2]]) == T)) > 0) {
    if(verbose)
      print("Filing missing values with multiple imputations:")
  
    tmp_ans <- mice(ans[,5:dim(ans)[2]], printFlag=ifelse(verbose, T, F))
    ans1 <- complete(tmp_ans)
    ans_final <- cbind(ans[,1:4], ans1)
  }
  
  if(verbose)
    print("Making final table...")
  ndim <- length(covariates)
  averages = matrix(nrow=1,ncol=length(covariates))
  
  dat <- ans_final[,1] #pid
  dat <- cbind(dat, ans_final[,2]) #sta (outcome)
  dat <- cbind(dat, ans_final[,3]) #tt1 (t1)
  dat <- cbind(dat, ans_final[,4]) #tt3 (t2)
  
  
  j <- 0
  i <- 0
  for(i in 0:(length(covariates)-1)) {
    dat <- cbind(dat, ans_final[,(5+i)]) 
    dat[2:dim(dat)[1],(5+j)] <- dat[1:(dim(dat)[1]-1),(5+j)]
    dat <- cbind(dat, ans_final[,(5+i)]) 
    averages[1,(i+1)] = dat[1,(5+j)]
    j <- j + 2
  }
  # Database should be in appropriate format:
  pid=dat[1,1]
  starttime = c(dat[1,3])
  for(i in 1:dim(dat)[1]) {
    if(dat[i,1] != pid) {
      avg <- c()
      for(ii in 0:(ndim-1)) {
        dat[(i+2),(5+ii)] = dat[(i+1),(6+ii)]
        avg <- c(avg, dat[i,(5+ii)])
        ii <- ii + 2
      }
      averages <- rbind(averages,avg)
      pid = dat[i,1]
      starttime <- c(starttime, dat[i,3])
    }
    if(dat[i,2] > 1) {
      dat[i,2] <- 1
    }
  }
  
  
  list(ans,dat)
}