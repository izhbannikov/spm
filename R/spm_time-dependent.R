# returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)

# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

optimize <- function(data, starting_params,  formulas, verbose, 
                     lower_bound, upper_bound, 
                     factr=1e-16, lmult=0.5,umult=2) {
  ### for test purposes:
  #data <- data1[[1]][,2:6]
  #starting_params <- list(a=-0.5, f1=80, Q=2e-8, f=80, b=5, mu0=1e-5)
  #formulas <- list(at="a1*t+a2", f1t="f1+f33*t", Qt="Q*exp(theta*t)", ft="f*t+ff", bt="b", mu0t="mu0*exp(theta*t)")
  #verbose=TRUE
  ### End of test purposes
  
  final_res <- list()
  # Current results:
  results <- list()
  
  N <- dim(data)[1]
  at <- NULL
  f1t <- NULL
  Qt <- NULL
  ft <- NULL
  bt <- NULL
  mu0t <- NULL
  variables <- c()
  
  # Assigning parameters:
  p.const.ind <- c()
  p.coeff.ind <- c()
  
  
  #---
  parameters <- trim(unlist(strsplit(formulas$at,"[\\+\\*\\(\\)]",fixed=F)))
  parameters <- parameters[which(!(parameters %in% c("t","exp")))]
  for(p in parameters) {assign(p,NULL, envir = .GlobalEnv); variables <- c(variables, p);}
  variables <- unique(variables)
  p.constants <- c()
  p.coeffs <- c()
  p.temp <- trim(unlist(strsplit(formulas$at,"[\\+\\-\\(\\)]",fixed=F)))
  for(i in 1:length(p.temp)) {
    if(grepl('t', p.temp[i]) == FALSE) {
      p.constants <- c(p.constants, p.temp[i])
    } else {
      p.temp.coeff <- trim(unlist(strsplit(p.temp[i],"[\\*]",fixed=F)))
      p.coeffs <- c(p.coeffs, p.temp.coeff[which(!(p.temp.coeff %in% c("t","exp")))])
    }
  }
  
  for(i in 1:length(p.constants)) {
    p.const.ind <- c(p.const.ind, which(variables == p.constants[i]))
  }
  for(i in 1:length(p.coeffs)) {
    p.coeff.ind <- c(p.coeff.ind, which(variables == p.coeffs[i]))
  }
  
  #----
  parameters <- trim(unlist(strsplit(formulas$f1t,"[\\+\\*\\(\\)]",fixed=F)))
  parameters <- parameters[which(!(parameters %in% c("t","exp")))]
  for(p in parameters) {assign(p,NULL, envir = .GlobalEnv); variables <- c(variables, p);} 
  variables <- unique(variables)
  p.constants <- c()
  p.coeffs <- c()
  p.temp <- trim(unlist(strsplit(formulas$f1t,"[\\+\\-\\(\\)]",fixed=F)))
  for(i in 1:length(p.temp)) {
    if(grepl('t', p.temp[i]) == FALSE) {
      p.constants <- c(p.constants, p.temp[i])
    } else {
      p.temp.coeff <- trim(unlist(strsplit(p.temp[i],"[\\*]",fixed=F)))
      p.coeffs <- c(p.coeffs, p.temp.coeff[which(!(p.temp.coeff %in% c("t","exp")))])
    }
  }
  
  for(i in 1:length(p.constants)) {
    p.const.ind <- c(p.const.ind, which(variables == p.constants[i]))
  }
  for(i in 1:length(p.coeffs)) {
    p.coeff.ind <- c(p.coeff.ind, which(variables == p.coeffs[i]))
  }
  
  #---
  parameters <- trim(unlist(strsplit(formulas$Qt,"[\\+\\*\\(\\)]",fixed=F)))
  parameters <- parameters[which(!(parameters %in% c("t","exp")))]
  for(p in parameters) {assign(p,NULL, envir = .GlobalEnv); variables <- c(variables, p);}
  variables <- unique(variables)
  p.constants <- c()
  p.coeffs <- c()
  p.temp <- trim(unlist(strsplit(formulas$Qt,"[\\+\\-\\(\\)]",fixed=FALSE)))
  for(i in 1:length(p.temp)) {
    if(grepl("*t", p.temp[i]) == FALSE) {
        p.temp.coeff <- trim(unlist(strsplit(p.temp[i],"[\\*]",fixed=FALSE)))
        p.constants <- c(p.constants,p.temp.coeff[which(!(p.temp.coeff %in% c("t","exp")))])
    } else {
        p.temp.coeff <- trim(unlist(strsplit(p.temp[i],"[\\*]",fixed=FALSE)))
        p.coeffs <- c(p.coeffs, p.temp.coeff[which(!(p.temp.coeff %in% c("t","exp")))])
    }
  }
  
  for(i in 1:length(p.constants)) {
    p.const.ind <- c(p.const.ind, which(variables == p.constants[i]))
  }
  for(i in 1:length(p.coeffs)) {
    p.coeff.ind <- c(p.coeff.ind, which(variables == p.coeffs[i]))
  }
  
  #---
  parameters <- trim(unlist(strsplit(formulas$ft,"[\\+\\*\\(\\)]",fixed=F)))
  parameters <- parameters[which(!(parameters %in% c("t","exp")))]
  for(p in parameters) {assign(p,NULL, envir = .GlobalEnv); variables <- c(variables, p);}
  variables <- unique(variables)
  p.constants <- c()
  p.coeffs <- c()
  p.temp <- trim(unlist(strsplit(formulas$ft,"[\\+\\-\\(\\)]",fixed=F)))
  for(i in 1:length(p.temp)) {
    if(grepl('t', p.temp[i]) == FALSE) {
      p.constants <- c(p.constants, p.temp[i])
    } else {
      p.temp.coeff <- trim(unlist(strsplit(p.temp[i],"[\\*]",fixed=F)))
      p.coeffs <- c(p.coeffs, p.temp.coeff[which(!(p.temp.coeff %in% c("t","exp")))])
    }
  }
  
  for(i in 1:length(p.constants)) {
    p.const.ind <- c(p.const.ind, which(variables == p.constants[i]))
  }
  for(i in 1:length(p.coeffs)) {
    p.coeff.ind <- c(p.coeff.ind, which(variables == p.coeffs[i]))
  }
  
  #---
  parameters <- trim(unlist(strsplit(formulas$bt,"[\\+\\*\\(\\)]",fixed=F)))
  parameters <- parameters[which(!(parameters %in% c("t","exp")))]
  for(p in parameters) {assign(p,NULL, envir = .GlobalEnv); variables <- c(variables, p);}
  variables <- unique(variables)
  p.constants <- c()
  p.coeffs <- c()
  p.temp <- trim(unlist(strsplit(formulas$bt,"[\\+\\-\\(\\)]",fixed=F)))
  for(i in 1:length(p.temp)) {
    if(grepl('t', p.temp[i]) == FALSE) {
      p.constants <- c(p.constants, p.temp[i])
    } else {
      p.temp.coeff <- trim(unlist(strsplit(p.temp[i],"[\\*]",fixed=F)))
      p.coeffs <- c(p.coeffs, p.temp.coeff[which(!(p.temp.coeff %in% c("t","exp")))])
    }
  }
  
  for(i in 1:length(p.constants)) {
    p.const.ind <- c(p.const.ind, which(variables == p.constants[i]))
  }
  for(i in 1:length(p.coeffs)) {
    p.coeff.ind <- c(p.coeff.ind, which(variables == p.coeffs[i]))
  }
  
  #---
  parameters <- trim(unlist(strsplit(formulas$mu0t,"[\\+\\*\\(\\)]",fixed=F)))
  parameters <- parameters[which(!(parameters %in% c("t","exp")))]
  for(p in parameters) {assign(p,NULL, envir = .GlobalEnv); variables <- c(variables, p);}
  variables <- unique(variables)
  p.constants <- c()
  p.coeffs <- c()
  p.temp <- trim(unlist(strsplit(formulas$mu0t,"[\\+\\-\\(\\)]",fixed=F)))
  for(i in 1:length(p.temp)) {
    if(grepl("*t", p.temp[i]) == FALSE) {
      p.temp.coeff <- trim(unlist(strsplit(p.temp[i],"[\\*]",fixed=FALSE)))
      p.constants <- c(p.constants,p.temp.coeff[which(!(p.temp.coeff %in% c("t","exp")))])
    } else {
      p.temp.coeff <- trim(unlist(strsplit(p.temp[i],"[\\*]",fixed=FALSE)))
      p.coeffs <- c(p.coeffs, p.temp.coeff[which(!(p.temp.coeff %in% c("t","exp")))])
    }
  }
  
  for(i in 1:length(p.constants)) {
    p.const.ind <- c(p.const.ind, which(variables == p.constants[i]))
  }
  for(i in 1:length(p.coeffs)) {
    p.coeff.ind <- c(p.coeff.ind, which(variables == p.coeffs[i]))
  }
  
  ########=============
  p.const.ind <- unique(p.const.ind)
  p.coeff.ind <- unique(p.coeff.ind)
  
  #variables <- unique(variables)
  stpar <- rep(0, length(variables))
  for(i in 1:length(stpar)) {
    stpar[p.const.ind[i]] <- unlist(starting_params, use.names = FALSE)[i]
  }
  names(stpar) <- variables
  
  for(p in names(stpar)) {
    results[[p]] <- stpar[p]
    assign(p, stpar[p], envir = globalenv())
  }
  
  # Lower and upper boundaries calculation:
  lower_bound <- c()
  for(i in 1:length(stpar)) {
    if(stpar[[i]] == 0) {stpar[[i]] = 1e-5}
    lower_bound <- c(lower_bound, ifelse(stpar[[i]] < 0, stpar[[i]] + lmult*stpar[[i]], stpar[[i]] - lmult*stpar[[i]]))
  }
  
  upper_bound <- c()
  for(i in 1:length(stpar)) {
    if(stpar[[i]] == 0) {stpar[[i]] = 1e-5}
    upper_bound <- c(upper_bound, ifelse(stpar[[i]] < 0, stpar[[i]] - umult*stpar[[i]], stpar[[i]] + umult*stpar[[i]]))
  }
  
  
  comp_func_params <- function(astring, f1string, qstring, fstring, bstring, mu0string) {
    at <<- eval(bquote(function(t) .(parse(text = astring)[[1]])))
    f1t <<- eval(bquote(function(t) .(parse(text = f1string)[[1]]))) 
    Qt <<- eval(bquote(function(t) .(parse(text = qstring)[[1]])))
    ft <<- eval(bquote(function(t) .(parse(text = fstring)[[1]])))
    bt <<- eval(bquote(function(t) .(parse(text = bstring)[[1]])))
    mu0t <<- eval(bquote(function(t) .(parse(text = mu0string)[[1]])))
  }
  
  iteration <- 0
  L.prev <- 0
  
  maxlik_t <- function(data, params) {
    stopflag <- FALSE
    
    if(verbose)
      cat("Iteration: ", iteration, "\n")
    
    for(p in names(stpar)) {
      assign(p, params[[p]], envir = globalenv())
      results[[p]] <<- params[[p]]
      if(verbose)
        cat(paste(p, get(p)), " ")
    }
    
    
    sigma_sq <- function(t1, t2) {
      # t2 = t_{j}, t1 = t_{j-1}
      ans <- bt(t1)*(t2-t1)
      ans
    }
  
    m <- function(y, t1, t2) {
      # y = y_{j-1}, t1 = t_{j-1}, t2 = t_{j}
      ans <- y + at(t1)*(y - f1t(t1))*(t2 - t1)
      ans
    }
  
    mu <- function(y, t) {
      ans <- mu0t(t) + (y - ft(t))^2*Qt(t)
    }
    
    #for(i in 1:length(results)) {
    #  if(length(intersect(results[[i]],c(lower_bound[i], upper_bound[i]))) >= 2) {
    #    cat("Parameter", names(results)[i], "achieved lower/upper bound. Process stopped.\n")
    #    cat(results[[i]],"\n")
    #    stopflag <- TRUE
    #    break
    #  }
    #}
  
    L <- 0
    
    if(stopflag == FALSE) {
      
      for(i in 1:N) {
        delta <- data[i,1]
        t1 <- data[i, 2]; t2 <- data[i, 3]
        ind <- ifelse(is.na(data[i, 5]), 0, 1)
        log_s <- -1*(mu(data[i, 4], t2-t1))
        if(ind == 0) {
          L <- L + (1 -delta)*(-1*log_s) + delta*(1-log_s)
        } else {
          L <- L + 0.5*pi*(-1*log(sigma_sq(t1, t2)) - (data[i,5] - m(data[i,4], t1, t2))^2/(2*sigma_sq(t1, t2)))
        }
      }
      
      assign("results", results, envir=.GlobalEnv)
      
      iteration <<- iteration + 1
      L.prev <<- L
      
      
      if(verbose) {
        cat("\n")
        cat(paste("L", L.prev), "\n")
      }
    
    } else {
      
      cat("Optimization stopped. Parametes achieved lower or upper bound.\nPerhaps you need more data or these returned parameters might be enough.\n")
      print("###########################################################")
      L <- NA
    
    }
    
    L
  }

  comp_func_params(formulas$at, formulas$f1t, formulas$Qt, formulas$ft, formulas$bt, formulas$mu0t)
  
  # Check if function parameters are consistent to those provided in starting list:
  if(verbose) {
    print(names(stpar))
    print(names(starting_params))
    print(variables)
  }
  
  #if( setequal(names(starting_params), variables) == FALSE) {
  #  stop("Provided set of function parameters is not equal to that one provided in starting list or vise-versa.")
  #}
  
  if(verbose == TRUE) {
    cat("Variables:\n")
    cat(variables,"\n")
    cat("Functions:\n")
    print(at)
    print(f1t)
    print(Qt)
    print(ft)
    print(bt)
    print(mu0t)
  }
  
  ## Optimization:
  tryCatch(optim(par = unlist(stpar), 
                  fn=maxlik_t, dat = as.matrix(data), 
                 control = list(fnscale=-1, trace=TRUE, maxit=10000), #ndeps=replicate(factr,n=length(lower_bound))), 
                  method="L-BFGS-B",#, 
                  lower=lower_bound,
                  upper=upper_bound),
           error=function(e) {print(e)}, 
           finally=NA)
  
  final_res <- list(results)
  final_res
}

#' spm_time_dep : a function that can handle time-dependant coefficients:
#' @param x : input data table.
#' @param start : a list of starting parameters, default: llist(a=-0.5, f1=80, Q=2e-8, f=80, b=5, mu0=1e-5, theta=0.08),
#' @param formulas : a list of formulas that define age (time) - dependency. Default: list(at="a1*t+a2", f1t="f1", Qt="Q*exp(theta*t)", ft="f", bt="b", mu0t="mu0*exp(theta*t)")
#' @return optimal coefficients
#' @examples
#' library(spm)
#' #Data preparation:
#' N <- 1000
#' data <- simdata_time_dep(N=2500)
#' opt.par <- spm_time_dep(data[,2:6], formulas=list(at="a", f1t="f1", Qt="Q*exp(theta*t)", ft="f", bt="b", mu0t="mu0*exp(theta*t)"), start=list(a=-0.5, f1=80, Q=2e-8, f=80, b=5, mu0=1e-5, theta=0.08))
#' opt.par
spm_time_dep <- function(x, 
                         start=list(a1=-0.5, a2=0.2, f1=80, Q=2e-8, f=80, b=5, mu0=1e-5),
                         formulas=list(at="a1*t+a2", f1t="f1", Qt="Q*exp(theta*t)", ft="f", bt="b", mu0t="mu0*exp(theta*t)"), 
                         verbose=TRUE,
                         lower_bound=NULL, upper_bound=NULL, 
                         factr=1e-16, lmult=0.5, umult=2) {
  
  data <- x
  
  # Optimization:
  #res = optimize(data, start, formulas, verbose, lower_bound, upper_bound, factr)
  res = optimize(data, start, formulas, verbose, factr=1e-16, lmult=lmult, umult=umult)
  invisible(res)
}