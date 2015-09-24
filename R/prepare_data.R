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

# Preparing data for stochastic process model
#'Output values include:
#'1). Database (data table), simulated for (slow) continuous optimization with arbitrary intervals between observations.
#'2). Database (data table), simulated for (quick) discrete optimization with fixed intervals between each observation.
#'@param longdat A table with longitude records.
#'@param vitstat A table with vital statistics (mortality).
#'@param interval A number of breaks between observations for discrete simulation. Default = 1 (no breaks).
#'@param col.status A name of column containing status variable (0/1 which indicate alive/dead). 
#'@param col.id A name of column containing patient ID. This ID should be the same in both longdat and vitstat tables.
#'@param col.age A name of age column.
#'@param col.age.event A name of event column.
#'@param covariates A list of covariates.
#'@param verbose A verbosing output indicator. Default=TRUE.
#'@return A list of two elements: first element contains a data table for continuous case, with arbitrary intervals between observations  and 
#'second element contains a data table for a discrete case (fixed intervals between observations).
#'@examples
#'library(spm)
#'#Reading longitude data:
#'longdat <- read.csv(system.file("data","longdat.csv",package="spm"))
#'# Prepare data for optimization:
#'vitstat <- read.csv(system.file("data","vitstat.csv",package="spm"))
#'# Remove unneeded NAs:
#'longdat.nonan <- longdat[which(is.na(longdat$Age) == F),]
#'vitstat.nonan <- vitstat[which(is.na(vitstat$BirthCohort) == F),]
#'data=prepare_data(longdat=longdat.nonan, vitstat=vitstat.nonan,interval=1, col.status="IsDead", col.id="ID", col.age="Age", col.age.event="LSmort", covariates=c("DBP"), verbose=T)
#'# Parameters estimation:
#'pars=spm(data,k = 1)
#'pars

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
  if((interval == 0) || (interval > 1)) {
    interval <- 1
  }
  
  
  #-----------Done parsing imput parameters---------------------#
  # Prepare data for continuous optimisation:
  data_cont <- prepare_data_cont(longdat, vitstat, col.status, col.id, col.age, col.age.event, covariates, verbose)
  
  # Prepare data for fast discrete optimization:
  data_discr <- prepare_data_discr(longdat, vitstat, interval, col.status, col.id, col.age, col.age.event, covariates, verbose)
  
  list(data_cont, data_discr)
}

prepare_data_cont <- function(longdat, vitstat, col.status, col.id, col.age, col.age.event, covariates, verbose) {
  #longdat=longdat.nonan
  #vitstat=vitstat.nonan
  #col.status="IsDead"
  #col.id="SubjID"
  #col.age="Age"
  #col.age.event="LSmort"
  #covariates=c("DBP", "BMI")
  #verbose=T
  
  # Split records by ID:
  prep.dat <- matrix(ncol=(4+2*length(covariates)),nrow=0)
  splitted <- split(longdat, longdat[[col.id]])
  vitstat.splitted <- split(vitstat, vitstat[[col.id]])
  
  for(iii in 1:length(splitted)) {
    nrows <- length(splitted[[iii]][[col.id]])
    id <- splitted[[iii]][[col.id]]
    case <- rep(0, nrows)
    case[nrows] <- vitstat.splitted[[iii]][[col.status]]
    t1 <- splitted[[iii]][[col.age]]
    t2 <- c(splitted[[iii]][[col.age]][-1], vitstat.splitted[[iii]][[col.age.event]])
    
    tmp.frame <- cbind(id, case, t1, t2)
    # Adding covariates:
    for(name in covariates) {
      tmp.frame <- cbind(tmp.frame, 
                         splitted[[iii]][[name]], 
                         c(splitted[[iii]][[name]][-1], NA))
      
    }
    prep.dat <- rbind(prep.dat, tmp.frame)
  }
  
  
  prep.dat <- prep.dat[rowSums( matrix(is.na(prep.dat[,5:dim(prep.dat)[2]]), ncol=2*length(covariates),byrow=T)) !=2*length(covariates),]
  prep.dat <- prep.dat[which(is.na(prep.dat[,4])==F),]
  head(prep.dat)  
  ans_final <- prep.dat
  if(length(which(is.na(prep.dat[,5:dim(prep.dat)[2]]) == T)) > 0) {
    if(verbose)
      cat("Filing missing values with multiple imputations:\n")
    
    tmp_ans <- mice(prep.dat[,5:dim(prep.dat)[2]], printFlag=ifelse(verbose, T, F),m = 2, maxit=2)
    ans1 <- complete(tmp_ans)
    #ans_final <- cbind(prep.dat[,1:3], ans1)
    ans_final <- cbind(prep.dat[,1:4], ans1)
  }
  
  if(verbose)
    cat("Making final table...\n")
  
  # Database should be in appropriate format:
  for(i in 1:(dim(ans_final)[1])) {
    if(ans_final[i,2] > 1) {
      ans_final[i,2] <- 1
    }
  }
  
  # Finalizing:
  ans_final <- ans_final[which(ans_final[,3] != ans_final[,4]),] # t1 must be different from t3
  ans_final <- ans_final[which(ans_final[,3] < ans_final[,4]),] # t1 must be less than t3
  # t1 must be equal t3 on previous step, if status = 0 and id is the same
  for(i in 2:dim(ans_final)[1]) {
    if((ans_final[i,3] != ans_final[(i-1),4]) & (ans_final[i,2] == 0) & (ans_final[i,1] == ans_final[(i-1),1])) {
      ans_final[i,3] <- ans_final[(i-1),4]
    }
  }
  
  colnames(ans_final) <- c("id", "case", "t1", "t2", unlist(lapply(1:length(covariates), function(n) {c(covariates[n], 
                                                                                                  paste(covariates[n],".next",sep=""))} )) )
  ans_final
  
}

prepare_data_discr <- function(longdat, vitstat, interval, col.status, col.id, col.age, col.age.event, covariates, verbose) {
  #longdat=longdat.nonan
  #vitstat=vitstat.nonan
  #interval = 3
  #col.status="IsDead"
  #col.id="ID"
  #col.age="Age"
  #col.age.event="LSmort"
  #covariates=c("DBP", "BMI")
  #verbose=T
  
  # Interpolation
  dt <- interval
  tt <- matrix(nrow=0, ncol=4)
  par <- matrix(nrow=0, ncol=length(covariates))
  
  # Split records by ID:
  splitted <- split(longdat, longdat[[col.id]])
  vitstat.splitted <- split(vitstat, vitstat[[col.id]])
  
  # For each particular person's record:
  for(iii in 1:length(splitted)) {
    if(!is.na(vitstat.splitted[[iii]][[col.age.event]]) & !is.na(vitstat.splitted[[iii]][[col.status]]) ) {
      
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
  
  ans <- ans[rowSums( matrix(is.na(ans[,5:dim(ans)[2]]), ncol=length(covariates),byrow=T)) !=length(covariates),]
  
  ans_final <- ans
  if(length(which(is.na(ans[,5:dim(ans)[2]]) == T)) > 0) {
    if(verbose)
      cat("Filing missing values with multiple imputations:\n")
    
    tmp_ans <- mice(ans[,5:dim(ans)[2]], printFlag=ifelse(verbose, T, F), m = 2, maxit = 2)
    ans1 <- complete(tmp_ans)
    ans_final <- cbind(ans[,1:4], ans1)
  }
  
  if(verbose)
    cat("Making final table...\n")
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
  for(i in 1:(dim(dat)[1]-1)) {
    if(dat[i,1] != pid) {
      for(ii in seq(0,(ndim-1),2)) {
        dat[(i+1),(5+ii)] = dat[i,(6+ii)]
      }
      pid = dat[i,1]
    }
    if(dat[i,2] > 1) {
      dat[i,2] <- 1
    }
  }
  
  colnames(dat) <- c("id", "case", "t1", "t2", unlist(lapply(1:length(covariates), function(n) {c(covariates[n], 
                                                                                                  paste(covariates[n],".next",sep="")
                                                                                                  )} 
                                                             )
                                                      ) 
                     )
  rownames(dat) <- 1:dim(dat)[1]
  dat
}

