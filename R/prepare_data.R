#'Data pre-processing for analysis with stochastic process model methodology.
#'@param x A path to the file with table of follow-up oservations (longitudinal table). 
#'File formats: csv, sas7bdat
#'@param col.id A name of column containing subject ID. 
#'This ID should be the same in both x (longitudinal) and y (vital statistics) tables.
#'None: if col.id not provided, the first column of the x and 
#'first column of the y will be used by default.
#'@param col.status A name of the column containing status variable 
#'(0/1, which is an indicator of death/censoring). 
#'Note: if not provided - then the column #2 from the y (vital statistics) dataset will be used.
#'@param col.age A name of age column (also called 't1'). 
#'This column represents a time (age) of measurement.
#'If not provided then the 3rd column from the longitudinal dataset (x) will be used.
#'@param col.age.event A name of 'event' column.
#'The event column indicates a time when the even occured (e.g. system failure).
#'Note: if not provided then the 3rd column from the y (vital statistics) dataset will be used.
#'@param covariates A list of covariates (physiological variables). 
#'If covariates not provided, then all columns from longitudinal table having index > 3 will be used as covariates. 
#'@param interval A number of breaks between observations for data for discrete model. 
#'This interval must be numeric (integer).
#'Default = 1 unit of time.
#'@param impute Multiple imputation ndicator. If TRUE then missing observations will be imputed with multiple imputation.
#'Default = TRUE.
#'@param verbose A verbosing output indicator. Default=FALSE.
#'@return A list of two elements: first element contains a preprocessed data for continuous model, with arbitrary intervals between observations  and 
#'second element contains a prepocessed data table for a discrete model (with constant intervals between observations).
#'@export
#'@examples \dontrun{ 
#'library(stpm) 
#'data <- prepare_data(x=system.file("data","longdat.csv",package="stpm"))
#'head(data[[1]])
#'head(data[[2]])
#'}
prepare_data <- function(x,
                         col.id=NULL, 
                         col.status=NULL,
                         col.age=NULL, 
                         col.age.event=NULL, 
                         covariates=NULL, 
                         interval=1, 
                         impute=TRUE,
                         verbose=FALSE) {
  
  
  if(interval < 1) {
    stop("Interval must be more or equal to 1.")
  } else if(interval != round(interval)) {
    stop("Interval must be integer.")
  }
  
  if(file_ext(x) == "csv") {
    #longdat <- read.csv(x)
    merged.data <- read.csv(x)
  } else if(file_ext(x) == "sas7bdat") {
    #longdat <- read.sas7bdat(x)
    merged.data <- read.sas7bdat(x)
  } else {
    stop(paste(x, ":", "unknown file format, it must be csv or sas7bdat."))
  }
  
  
  #if(file_ext(y) == "csv") {
  #  vitstat <- read.csv(y)
  #} else if(file_ext(y) == "sas7bdat") {
  #  vitstat <- read.sas7bdat(y)
  #} else {
  #  stop(paste(y, ":", "unknown file format, it must be csv or sas7bdat."))
  #}
  
  # Parsing input parameters in order to check for errors:
  if( !is.null(col.status) ) {
    if( !(col.status %in% colnames(merged.data)) ) {
      stop(paste("Status column",col.status, "not found in data table. Aborting."))
    }
  }
  
  if( !is.null(col.id) ) { 
    if( !(col.id %in% colnames(merged.data)) ) {
      stop(paste("ID column",col.id, "not found in data table. Aborting."))
    }
  }
  
  if( !is.null(col.age) ) {
    if( !(col.age %in% colnames(merged.data)) ) {
      stop(paste("Age column",col.age, "not found in longdat table. Aborting."))
    }
  }
  
  if( !is.null(col.age.event) ) { 
    if( !(col.age.event %in% colnames(merged.data)) ) {
      stop(paste("Event column",col.age.event, "not found in data table. Aborting."))
    }
  }
  
  if(!is.null(covariates)) {
    for(c in covariates) {
      if( !(c %in% colnames(merged.data)) ) {
        stop(paste("Covariate",c, "not found. Aborting."))
      }
    }
  } else if(is.null(covariates)) {
    col.covar.ind <- 4:dim(merged.data)[2]
  }
  
  if((interval == 0) || (interval < 1)) {
    interval <- 1
  }
  
  #-----------Done parsing imput parameters---------------------#
  
  #merged.data <- merge(x = longdat, y = vitstat, by.x = col.id, by.y=col.id)
  
  if(!is.null(col.status)) {
    col.status.ind <- grep(paste("\\b", col.status, "\\b", sep=""), colnames(merged.data))
  } else {
    col.status.ind <- 2
  }
  
  if(!is.null(col.id)) {
    col.id.ind <- grep(paste("\\b", col.id, "\\b", sep=""), colnames(merged.data))
  } else {
    col.id.ind <- 1
  }
  
  if(!is.null(col.age)) {
    col.age.ind <- grep(paste("\\b", col.age, "\\b", sep=""), colnames(merged.data))
  } else {
    col.age.ind <- 3
  }
  
  if(!is.null(col.age.event)) {
    col.age.event.ind <- grep(paste("\\b", col.age.event, "\\b", sep=""), colnames(merged.data))
  } else {
    #col.age.event.ind <- 3
    col.age.event.ind <- col.age.ind
  }
  
  if(!is.null(covariates)) {
    col.covar.ind <- c()
    for(c in covariates) {
        col.covar.ind <- c(col.covar.ind, grep(paste("\\b", c, "\\b", sep=""), colnames(merged.data)))
    }
  } else {
    #col.covar.ind <- 4:dim(longdat)[2]
    col.covar.ind <- 4:dim(merged.data)[2]
  }
  
  merged.data <- merged.data[which(!is.na(merged.data[ , col.age.ind])),]
  
  # Prepare data for continuous optimisation:
  data_cont <- prepare_data_cont(merged.data, col.status.ind, col.id.ind, col.age.ind, col.age.event.ind, col.covar.ind, verbose, impute, interval)
  
  # Prepare data for fast discrete optimization:
  data_discr <- prepare_data_discr(merged.data, interval, col.status.ind, col.id.ind, col.age.ind, col.age.event.ind, col.covar.ind, verbose, impute)
  
  list(model.continuous=data_cont, model.discrete=data_discr)
}

#'Prepares continuouts-time dataset.
#'@param merged.data a longitudinal study dataset.
#'@param col.status.ind index of "status" column.
#'@param col.id.ind subject id column index.
#'@param col.age.ind index of the age column.
#'@param col.age.event.ind an index of the column which represents the time in which event occured.
#'@param col.covar.ind a set of column indexes which represent covariates.
#'@param verbose turns on/off verbosing output.
#'@param impute Multiple imputation ndicator. If TRUE then missing observations will be imputed with multiple imputation.
#'@param dt interval between observations.
prepare_data_cont <- function(merged.data, 
                              col.status.ind, 
                              col.id.ind, 
                              col.age.ind, 
                              col.age.event.ind, 
                              col.covar.ind, 
                              verbose,
                              impute,
                              dt) {
  
  # Split records by ID:
  prep.dat <- data.frame(matrix(ncol=(4+2*length(col.covar.ind)),nrow=0))
  splitted <- split(merged.data, merged.data[ , col.id.ind])
  
  for(iii in 1:length(splitted)) {
    nrows <- length(splitted[[iii]][ , col.id.ind])
    id <- splitted[[iii]][ , col.id.ind]
    case <- splitted[[iii]][, col.status.ind]
    t1 <- splitted[[iii]][ , col.age.ind]
    t2 <- c(splitted[[iii]][ , col.age.ind][-1], tail(splitted[[iii]][ , col.age.event.ind],n=1))
    
    tmp.frame <- cbind(id, case, t1, t2)
    # Adding covariates:
    for(ind in col.covar.ind) {
      tmp.frame <- cbind(tmp.frame, 
                         splitted[[iii]][, ind], 
                         c(splitted[[iii]][, ind][-1], NA))
      
    }
    prep.dat <- rbind(prep.dat, tmp.frame)
    
  }
  
  
  #prep.dat <- prep.dat[rowSums( matrix(is.na(prep.dat[,5:dim(prep.dat)[2]]), ncol=2*length(covariates),byrow=T)) !=2*length(covariates),]
  prep.dat <- prep.dat[which(!is.na(prep.dat[,4])),]
  
  if(verbose) {
    head(prep.dat)  
  }
  
  ans_final <- prep.dat
  
  if(impute) {
    if(verbose)
      cat("Filing missing values with multiple imputations:\n")
    
    tmp_ans <- mice(prep.dat[,5:dim(prep.dat)[2]], printFlag=ifelse(verbose, TRUE, FALSE),m = 2, maxit=2)
    ans1 <- complete(tmp_ans)
    ans_final <- cbind(prep.dat[,1:4], ans1)
  }
  
  if(verbose)
    cat("Making final table...\n")
  
  ans_final <- ans_final[which(!is.na(ans_final$case)), ]
  
  # Database should be in appropriate format:
  for(i in 1:(dim(ans_final)[1])) {
    if(ans_final[i,2] > 1) {
      ans_final[i,2] <- 1
    }
  }
  
  # Finalizing:
  
  ans_final <- ans_final[which(ans_final[,3] <= ans_final[,4]),] # t1 must be less than t3
  # t1 must be equal t3 on previous step, if status = 0 and id is the same
  ndim <- dim(ans_final)[2] - 4
  for(i in 2:(dim(ans_final)[1]-1)) {
    if( (ans_final$case[i] == 0) & (ans_final$id[i] == ans_final$id[(i-1)]) & (ans_final$id[i] == ans_final$id[(i+1)]) ) {
      for(ii in seq(0,(ndim-1),2)) {
        ans_final[(i+1),(5+ii)] <- ans_final[i,(6+ii)]
      }
    } else if(ans_final$case[i] == 1) {
      
      for(ii in seq(0,(ndim-1),2)) {
        ans_final[i,(6+ii)] <- NA
      }
    } else if( ans_final$id[i] != ans_final$id[(i+1)] ) {
      if(ans_final$t1[i] >= ans_final$t2[i]) {
        ans_final$t2[i] <- ans_final$t1[i] + dt/2
      }
    }
    
    
  }
  
  colnames(ans_final) <- c("id", "case", "t1", "t2", unlist(lapply(1:length(col.covar.ind), function(n) {c(names(merged.data)[col.covar.ind[n]], 
                                                                                                           paste(names(merged.data)[col.covar.ind[n]],".next",sep=""))} )) )
  data.frame(ans_final)
  
}

#'Prepares discrete-time dataset.
#'@param merged.data a longitudinal study dataset.
#'@param interval interval between observations.
#'@param col.status.ind index of "status" column.
#'@param col.id.ind subject id column index.
#'@param col.age.ind index of the age column.
#'@param col.age.event.ind an index of the column which represents the time in which event occured.
#'@param col.covar.ind a set of column indexes which represent covariates.
#'@param verbose turns on/off verbosing output.
#'@param impute Multiple imputation ndicator. If TRUE then missing observations will be imputed with multiple imputation.
prepare_data_discr <- function(merged.data, interval, col.status.ind, col.id.ind, col.age.ind, col.age.event.ind, col.covar.ind, verbose, impute) {
  #---DEBUG---#
  #longdat <- read.sas7bdat("/Volumes/G/spm/data/covar_aric_gru.sas7bdat")
  #vitstat <- read.sas7bdat("/Volumes/G/spm/data/mortality_aric_all_gru.sas7bdat")
  #interval <- 1
  #col.id="SubjID"
  #col.age="Age"
  #col.age.event="LSmort"
  #col.status="IsDead"
  #covariates="BMI"
  #col.status.ind <- grep(paste("\\b", col.status, "\\b", sep=""), colnames(vitstat))
  #col.id.ind <- grep(paste("\\b", col.id, "\\b", sep=""), colnames(vitstat)) 
  #col.age.ind <- grep(paste("\\b", col.age, "\\b", sep=""), colnames(longdat))
  #col.age.event.ind <- grep(paste("\\b", col.age.event, "\\b", sep=""), colnames(vitstat))
  #col.covar.ind <-  grep(paste("\\b", "BMI", "\\b", sep=""), colnames(longdat))
  #verbose <- TRUE
  #longdat <- longdat[which(!is.na(longdat[ , col.age.ind])),]
  
  
  #longdat = longdat
  #vitstat = vitstat
  #col.id = "ID"
  #col.status = "IsDead"
  #col.age = "Age"
  #col.age.event = "LSmort"
  #covariates = "DBP"
  #
  #col.status.ind = col.status.ind
  #col.id.ind  = col.id.ind
  #col.age.ind = col.age.ind
  #col.age.event.ind = col.age.event.ind
  #col.covar.ind = col.covar.ind
  #---END DEBUG---#
  
  #'Filling the last cell
  fill_last <- function(x) {
    na_idx <- which(is.na(x))
    unique_elements <- unique(x[-na_idx])
    set_diff <- unique_elements[length(unique_elements)]
    x[na_idx] <- set_diff
    x
  }
  
  
  # Interpolation
  dt <- interval
  tt <- matrix(nrow=0, ncol=4)
  par <- matrix(nrow=0, ncol=length(col.covar.ind))
  
  # Split records by ID:
  splitted <- split(merged.data, merged.data[, col.id.ind])
  
  # For each particular person's record:
  for(iii in 1:length(splitted)) {
    
    if( !is.na(tail(splitted[[iii]][ , col.age.event.ind],n=1)) & !is.na(tail(splitted[[iii]][ , col.status.ind],n=1)) ) {
      if(verbose) {
        print(paste(iii, "individual processed."))
      }
      # Individual ID:
      id <- splitted[[iii]][ , col.id.ind][1]
      nrows <- (tail(splitted[[iii]][ , col.age.ind], n=1) - floor(splitted[[iii]][ , col.age.ind][1]))/dt + 1
      
      # Perform approximation using two points:
      t1.approx <- matrix(ncol=4, nrow=nrows)
      t1.approx[,1] <- id
      t1.approx[,2] <- 0
      t1.approx[nrows,2] <- tail(splitted[[iii]][ , col.status.ind],n=1) #Last value
      t1.approx[,3] <- seq(floor(splitted[[iii]][ , col.age.ind][1]), splitted[[iii]][ , col.age.ind][length(splitted[[iii]][ , col.age.ind])], by=dt)
      if(nrows > 1) {
        t1.approx[,4] <- c(t1.approx[,3][2:nrows], tail(splitted[[iii]][ , col.age.event.ind],n=1))
      } else {
        t1.approx[,4] <- tail(splitted[[iii]][ , col.age.event.ind],n=1)
      }
      
      tt <- rbind(tt,t1.approx)
      par1.approx <- matrix(ncol=length(col.covar.ind), nrow=nrows, NA)
      
      j <- 1
      for(ind in col.covar.ind) {
        if ( (length(splitted[[iii]][, ind]) > 1) & (length(which(!is.na(splitted[[iii]][, ind]))) > 0) ) {
          if(length(which(!is.na(splitted[[iii]][, ind]))) == 1) {
            splitted[[iii]][, ind] <- fill_last(splitted[[iii]][, ind])
          }
          # Fill NAs by linear approximation with approx():
          #nn <- length(splitted[[iii]][, ind])
          #splitted[[iii]][, ind] <- approx(splitted[[iii]][, ind],n=nn)$y
          
          #aprx <- c()
          #for(k in 1:(length(splitted[[iii]][, col.age.ind])-1)) {
          #  nr <- ceiling((splitted[[iii]][, col.age.ind][k+1] - splitted[[iii]][, col.age.ind][k])/dt) + 1
          #  aprx <- c(aprx[1:length(aprx)-1], approx(splitted[[iii]][, ind][k:(k+1)], n=nr)$y)
          #  print(nr)
          #  print(aprx)
          #}
          aprx <- approx(splitted[[iii]][, ind],n=nrows)$y
          #print(par1.approx[,j])
          par1.approx[,j] <- aprx
          #par1.approx[,j] <-  approx(splitted[[iii]][, ind], n=nrows)$y
        }
        
        j <- j + 1
        
      }
      par <- rbind(par,par1.approx)
    }
  }
  
  ans <- cbind(tt,par)
  colnames(ans) <- c("id", "case", "t1", "t2", names(merged.data)[col.covar.ind])
  
  ans <- data.frame(ans[rowSums( matrix(is.na(ans[,5:dim(ans)[2]]), ncol=length(col.covar.ind),byrow=T)) !=length(col.covar.ind),])
  
  ans_final <- ans
  
  if(verbose)
    cat("Making final table...\n")
  
  ndim <- length(col.covar.ind)
  averages = matrix(nrow=1,ncol=length(col.covar.ind))
  
  dat <- ans_final[,1] #pid
  dat <- cbind(dat, ans_final[,2]) #sta (outcome)
  dat <- cbind(dat, ans_final[,3]) #tt1 (t1)
  dat <- cbind(dat, ans_final[,4]) #tt3 (t2)
  
  
  for(i in 0:(length(col.covar.ind)-1)) {
    dat <- cbind(dat, ans_final[,(5+i)]) 
    dat <- cbind(dat, ans_final[,(5+i)])
  }
  
  if(impute) {
    if(verbose)
      cat("Filing missing values with multiple imputations:\n")
    
    tryCatch({
      tmp_ans <- mice(as.data.frame(dat[,5:dim(dat)[2]]), printFlag=ifelse(verbose, TRUE, FALSE), m = 2, maxit = 2)
      ans1 <- complete(tmp_ans)
      ans_final <- cbind(ans[,1:4], ans1)
      dat <- ans_final
    }, error=function(e) {
      print(e)
      
    })
  }
  
  
  k <- 0
  for(i in 0:(length(col.covar.ind)-1)) {
    for(j in 1:(dim(dat)[1]-1)) {
      if(dat[j,1] != dat[(j+1), 1]) {
        dat[j, (5+k+1)] <- NA
        
        if(dat[j,3] >= dat[j,4]) {
          dat[j,4] <- dat[j,3] + dt/2
        }
        
      } else {
        dat[j, (5+k+1)] <- dat[(j+1), (5+k)]
      }
      
      if(dat[j,2] > 1) {
        dat[j,2] <- 1
      }
      
    }
    k <- k+2
  }
  
  colnames(dat) <- c("id", "case", "t1", "t2", 
                     unlist(lapply(1:length(col.covar.ind), 
                                   function(n) {c(names(merged.data)[col.covar.ind[n]], 
                                                paste(names(merged.data)[col.covar.ind[n]],".next",sep=""))})))
  rownames(dat) <- 1:dim(dat)[1]
  return(data.frame(dat))
}
