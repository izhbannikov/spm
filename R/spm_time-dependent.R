#'Returns string w/o leading whitespace
#'@param x a string to trim
trim.leading <- function (x)  sub("^\\s+", "", x)

#'Returns string w/o trailing whitespace
#'@param x a string to trim
trim.trailing <- function (x) sub("\\s+$", "", x)

#'Returns string w/o leading or trailing whitespace
#'@param x a string to trim
trim <- function (x) gsub("^\\s+|\\s+$", "", x)


#'spm_time_dep : a function for the model with 
#'time-dependent model parameters.
#'@param x Input data table.
#'@param start A list of starting parameters, default: 
#'\code{start=list(a=-0.5, f1=80, Q=2e-8, f=80, b=5, mu0=1e-5)}.
#'@param frm A list of formulas that define age (time) - dependency. 
#'Default: \code{frm=list(at="a", f1t="f1", Qt="Q", ft="f", bt="b", 
#'mu0t="mu0")}.
#'@param stopifbound Estimation stops if at least one parameter 
#'achieves lower or upper boundaries. Default: \code{FALSE}.
#'@param lb Lower bound of parameters under estimation.
#'@param ub Upper bound of parameters under estimation.
#'@param verbose Turns on verbosing output.
#'@param opts A list of options for \code{nloptr}.
#'Default value: \code{opt=list(algorithm="NLOPT_LN_NELDERMEAD", 
#'maxeval=100, ftol_rel=1e-8)}.
#'Please see \code{nloptr} documentation for more information.
#'@return A set of estimates of \code{a}, \code{f1}, \code{Q}, 
#'\code{f}, \code{b}, \code{mu0}.
#'@return status Optimization status (see documentation for nloptr package).
#'@return LogLik A logarithm likelihood.
#'@return objective A value of objective function (given by nloptr).
#'@return message A message given by \code{nloptr} optimization function 
#'(see documentation for nloptr package).
#'@references Yashin, A. et al (2007), Health decline, 
#'aging and mortality: how are they related? 
#'Biogerontology, 8(3), 291-302.<DOI:10.1007/s10522-006-9073-3>.
#'@export
#'@examples
#'library(stpm)
#'#Data preparation:
#'n <- 5
#'data <- simdata_time_dep(N=n)
#'# Estimation:
#'opt.par <- spm_time_dep(data)
#'opt.par
#'
spm_time_dep <- function(x, 
                         start=list(a=-0.05, f1=80, Q=2e-8, f=80, b=5, mu0=1e-3),
                         frm=list(at="a", f1t="f1", Qt="Q", ft="f", 
                                  bt="b", mu0t="mu0"), 
                         stopifbound=FALSE, 
                         lb=NULL, ub=NULL,
                         verbose=FALSE, opts=list(algorithm="NLOPT_LN_NELDERMEAD", 
                                                  maxeval=100, ftol_rel=1e-8)) {
  
  #avail_algorithms <- c("NLOPT_GN_DIRECT", "NLOPT_GN_DIRECT_L",
  #                      "NLOPT_GN_DIRECT_L_RAND", "NLOPT_GN_DIRECT_NOSCAL",
  #                      "NLOPT_GN_DIRECT_L_NOSCAL",
  #                      "NLOPT_GN_DIRECT_L_RAND_NOSCAL",
  #                      "NLOPT_GN_ORIG_DIRECT", "NLOPT_GN_ORIG_DIRECT_L",
  #                      "NLOPT_GN_CRS2_LM", 
  #                      "NLOPT_LN_COBYLA", "NLOPT_LN_NEWUOA",
  #                      "NLOPT_LN_NEWUOA_BOUND", "NLOPT_LN_NELDERMEAD",
  #                      "NLOPT_LN_SBPLX",
  #                      "NLOPT_LN_BOBYQA", "NLOPT_GN_ISRES")
  
  if(is.null(start)) {
      start <- list(a=-0.05, f1=80, Q=2e-8, f=80, b=5, mu0=1e-3)
  }
  
  #--------------Begin of optimize function-------------------#
  optimize <- function(data, starting_params,  formulas, verbose, 
                       lb, ub, 
                       stopifbound,
                       opts) {
    
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
    #for(p in parameters) {assign(p,NULL, envir = .GlobalEnv); variables <- c(variables, p);}
    for(p in parameters) {assign(p,NULL, envir=baseenv()); variables <- c(variables, p);}
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
    #for(p in parameters) {assign(p,NULL, envir = .GlobalEnv); variables <- c(variables, p);} 
    for(p in parameters) {assign(p,NULL, envir=baseenv()); variables <- c(variables, p);} 
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
    #for(p in parameters) {assign(p,NULL, envir = .GlobalEnv); variables <- c(variables, p);}
    for(p in parameters) {assign(p,NULL, envir=baseenv()); variables <- c(variables, p);}
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
    #for(p in parameters) {assign(p,NULL, envir = .GlobalEnv); variables <- c(variables, p);}
    for(p in parameters) {assign(p,NULL, envir=baseenv()); variables <- c(variables, p);}
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
    for(p in parameters) {assign(p,NULL, envir=baseenv()); variables <- c(variables, p);}
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
    for(p in parameters) {assign(p,NULL, envir=baseenv()); variables <- c(variables, p);}
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
    
    variables <- variables[which(variables != "0")]
    
    
    stpar <- rep(0, length(variables))
    for(i in 1:length(stpar)) {
      stpar[i] <- unlist(starting_params, use.names = FALSE)[i]
    }
    names(stpar) <- variables
    
    for(p in names(stpar)) {
      results[[p]] <- stpar[p]
      assign(p, stpar[p], envir=baseenv())
    }
    
    # Lower and upper boundaries calculation:
    
    lower_bound <- c()
    
    if(is.null(lb)) {
      for(i in 1:length(stpar)) {
        if(stpar[[i]] == 0) {stpar[[i]] = 1e-12}
        lower_bound <- c(lower_bound, ifelse(stpar[[i]] < 0, stpar[[i]] + 0.25*stpar[[i]], stpar[[i]] - 0.25*stpar[[i]]))
      }
    } else {
      lower_bound <- lb
    }
    
    upper_bound <- c()
    if(is.null(ub)) {
      for(i in 1:length(stpar)) {
        if(stpar[[i]] == 0) {stpar[[i]] = 1e-12}
        upper_bound <- c(upper_bound, ifelse(stpar[[i]] < 0, stpar[[i]] - 0.25*stpar[[i]], stpar[[i]] + 0.25*stpar[[i]]))
      }
    } else {
      upper_bound <- ub
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
    
    maxlik_t <- function(params) {
      stopflag <- FALSE
      
      if(verbose)
        cat("Iteration: ", iteration, "\n")
      
      names(params) <- names(stpar)
      
      for(p in names(stpar)) {
        assign(p, params[[p]], envir=baseenv())
        results[[p]] <<- params[[p]]
        if(verbose)
          cat(paste(p, results[[p]]), " ")
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
      
      mu <- function(y,t) {
        ans <- mu0t(t) + (y - ft(t))^2*Qt(t)
        ans
      }
      
      if(stopifbound) {
        for(i in 1:length(results)) {
          if(length(intersect(results[[i]],c(lower_bound[i], upper_bound[i]))) >= 2) {
            cat("Parameter", names(results)[i], "achieved lower/upper bound. Process stopped.\n")
            cat(results[[i]],"\n")
            stopflag <- TRUE
            break
          }
        }
      }
      
      L <- 0
      
      if(stopflag == FALSE) {
        
        for(i in 1:N) {
          delta <- data[i,1]
          t1 <- data[i, 2]; t2 <- data[i, 3]
          if(t1 < t2) {
            ind <- ifelse(is.na(data[i, 5]), 0, 1)
            S <- exp(-1*mu(data[i, 4],t1)*(t2-t1))
            if(S <= 1e-5) {
              S <- 1e-5
            }
            if(ind == 0) {
              L <- L + (1 - delta)*log(S) + delta*log(1-S)
              #if(is.nan(L)) {
              #  print(mu(data[i, 4],t1)*(t2-t1))
              #  print(S)
              #  print(data[i,])
              #  stop()
              #}
            } else {
              yj <- data[i,5]
              mj <- m(data[i,4], t1, t2)
              pn <- -1*log(sqrt(2*pi)*sqrt(sigma_sq(t1, t2))) - (yj - mj)^2/(2*sigma_sq(t1, t2))
              L <- L + pn + (1 - delta)*log(S) + delta*log(1-S)
            }
            
          }
        }
        
        assign("results", results, envir=baseenv())
        
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
      
      L <- -1*L
      L
    }
    
    comp_func_params(formulas$at, formulas$f1t, formulas$Qt, formulas$ft, formulas$bt, formulas$mu0t)
    
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
    
    # Optimization:
    if(verbose) {
      cat("Lower bound:\n")
      print(lower_bound)
      cat("Upper bound:\n")
      print(upper_bound)
      
    }
    
    
    tryCatch({ans <- nloptr(x0 = unlist(stpar), 
                            eval_f = maxlik_t, opts = opts,
                            lb = lower_bound, ub = upper_bound)
    i <- 1
    for(p in names(stpar)) {
      results[[p]] <<- ans$solution[i]
      i <- i + 1
      if(verbose)
        cat(paste(p, get(p)), " ")
    }
    results[["status"]] <- ans$status
    results[["LogLik"]] <- L.prev
    results[["objective"]] <- ans$objective
    results[["message"]] <- ans$message
    },  
    error=function(e) {if(verbose  == TRUE) {print(e)}}, 
    finally=NA)
    
    final_res <- list(results)
    final_res
  }
  #---------------------End of optimize---------------------------#
  
  formulas <- frm
  data <- as.matrix(x)
  data <- data[, 2:dim(data)[2]]
  formulas.work = list(at="a", f1t="f1", Qt="Q", ft="f", bt="b", mu0t="mu0")
  for(item in names(formulas)) {
    formulas.work[[item]] <- formulas[[item]]
  }
  # Optimization:
  res = optimize(data, start, formulas.work, verbose, lb, ub, stopifbound, opts)
  invisible(res)
}