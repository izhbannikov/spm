# Preparing data for stochastic process model
prepare_data <- function(interval=1, col.status="IndicatorDeath", col.id="ID", col.age="Age", col.age.next="AgeNext", covariates=c("DBP", "BMI", "DBP1", "DBP2", "Weight", "Height") ) {
  
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
  
  splitted <- split(longdat.nonan, longdat.nonan[[col.id]])
  
  for(iii in 1:length(splitted)) {
    print(iii)
    id <- splitted[[iii]][[col.id]][1]
    nrows <- (tail(splitted[[iii]][[col.age]], n=1) - splitted[[iii]][[col.age]][1])/dt + 1
    # Perform approximation:
    t1.approx <- matrix(ncol=4, nrow=nrows)
    t1.approx[,1] <- id
    t1.approx[,2] <- 0
    t1.approx[nrows,2] <- splitted[[iii]][[col.status]][length(splitted[[iii]][[col.status]])] #Last value
    t1.approx[,3] <- seq(splitted[[iii]][[col.age]][1], splitted[[iii]][[col.age]][length(splitted[[iii]][[col.age]])], by=dt)
    if(nrows > 1) {
      t1.approx[,4] <- c(t1.approx[,3][2:nrows], splitted[[iii]][[col.age.next]][length(splitted[[iii]][[col.age.next]])])
    } else {
      t1.approx[,4] <- splitted[[iii]][[col.age.next]][1]
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
  
  ans=cbind(tt,par)
  colnames(ans) <- c("ID", "CASE", "T1", "T3", covariates)
  ans <- ans[rowSums(is.na(ans[,5:dim(ans)[2]]))!=length(covariates),]
  #ans <- ans[which(is.na(ans)==F,arr.ind=TRUE)[,1],]
  
  print("Filing missing values with multiple imputations:")
  tmp_ans=mice(ans)
  ans_final = complete(tmp_ans)
  ans_final
}
