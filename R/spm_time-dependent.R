# returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)

# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

optimize <- function(data, starting_params,  formulas, verbose) {
  N <- dim(data)[1]
  at <- NULL#; a <- starting_params$a;
  f1t <- NULL#; f1 <- starting_params$f1;
  Qt <- NULL#; Q <- starting_params$Q;
  ft <- NULL#; f <- starting_params$f;
  bt <- NULL#; b <- starting_params$b;
  mu0t <- NULL#; mu0 <- starting_params$mu0;
  #theta <- starting_params$theta
  variables <- c()
  # Assigning parameters:
  #---
  parameters = trim(unlist(strsplit(formulas$at,"[\\+\\*\\(\\)]",fixed=F)))
  parameters = parameters[which(!(parameters %in% c("t","exp")))]
  for(p in parameters) {assign(p,NULL, envir = .GlobalEnv); variables <- c(variables, p);}
  
  #----
  parameters = trim(unlist(strsplit(formulas$f1t,"[\\+\\*\\(\\)]",fixed=F)))
  parameters = parameters[which(!(parameters %in% c("t","exp")))]
  for(p in parameters) {assign(p,NULL, envir = .GlobalEnv); variables <- c(variables, p);}  
  #---
  parameters = trim(unlist(strsplit(formulas$Qt,"[\\+\\*\\(\\)]",fixed=F)))
  parameters = parameters[which(!(parameters %in% c("t","exp")))]
  for(p in parameters) {assign(p,NULL, envir = .GlobalEnv); variables <- c(variables, p);}
  
  #---
  parameters = trim(unlist(strsplit(formulas$ft,"[\\+\\*\\(\\)]",fixed=F)))
  parameters = parameters[which(!(parameters %in% c("t","exp")))]
  for(p in parameters) {assign(p,NULL, envir = .GlobalEnv); variables <- c(variables, p);}
  #---
  parameters = trim(unlist(strsplit(formulas$bt,"[\\+\\*\\(\\)]",fixed=F)))
  parameters = parameters[which(!(parameters %in% c("t","exp")))]
  for(p in parameters) {assign(p,NULL, envir = .GlobalEnv); variables <- c(variables, p);}
  #---
  parameters = trim(unlist(strsplit(formulas$mu0t,"[\\+\\*\\(\\)]",fixed=F)))
  parameters = parameters[which(!(parameters %in% c("t","exp")))]
  for(p in parameters) {assign(p,NULL, envir = .GlobalEnv); variables <- c(variables, p);}
  
  variables <- unique(variables)
  
  for(p in names(starting_params)) {
    assign(p, starting_params[[p]], envir = globalenv())
  }
  #print(variables)
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
    
    if(verbose)
      cat("Iteration: ", iteration, "\n")
    
    for(p in names(starting_params)) {
      assign(p, params[[p]], envir = globalenv())
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
  
    L <- 0
  
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
    
    iteration <<- iteration + 1
    L.prev <<- L
    
    if(verbose) {
      cat("\n")
      cat(paste("L", L.prev), "\n")
    }
    
    L
  }

  comp_func_params(formulas$at, formulas$f1t, formulas$Qt, formulas$ft, formulas$bt, formulas$mu0t)
  
  ## Optimization:
  #result <- optim(par = unlist(starting_params), 
  #              fn=maxlik_t, dat = as.matrix(data), control = list(fnscale=-1, trace=T, maxit=10000, factr=1e-16, ndeps=c(1e-12, 1e-12, 1e-12,1e-12,1e-16,1e-12,1e-12,1e-12,1e-12)), 
  #              method="L-BFGS-B", lower = c(-0.5, -0.5, -1, 0,1e-12,1e-6,1e-6,1e-6, 1e-4), upper = c(0, 3, 0, Inf, 1e-7, Inf, Inf, 1, 0.1))
  print(unlist(starting_params))
  result <- optim(par = unlist(starting_params), 
                  fn=maxlik_t, dat = as.matrix(data), control = list(fnscale=-1, trace=T, maxit=10000, factr=1e-16), 
                  method="L-BFGS-B")
  
  result

}

#' spm_time_dep : a function that can handle time-dependant coefficients:
#' @param start : a list of starting parameters, default: llist(a=-0.5, f1=80, Q=2e-8, f=80, b=5, mu0=1e-5, theta=0.08),
#' @param formulas : a list of formulas that define age (time) - dependency. Default: list(at="a1*t+a2", f1t="f1", Qt="Q*exp(theta*t)", ft="f", bt="b", mu0t="mu0*exp(theta*t)")
#' @return optimal coefficients
#' @examples
#' library(spm)
#' Data preparation:
#' N <- 1000
#' data <- simdata_cont(N=N, aH=-0.05, f1H=80, QH=2e-8, fH=80, bH=5, mu0H=2e-5, thetaH=0.08)
#' opt.par <- spm_time_dep(data[,2:6], formulas=list(at="a", f1t="f1", Qt="Q*exp(theta*t)", ft="f", bt="b", mu0t="mu0*exp(theta*t)"), start=list(a=-0.5, f1=80, Q=2e-8, f=80, b=5, mu0=1e-5, theta=0.08))
#' opt.par
spm_time_dep <- function(data, 
                         start=list(a1=-0.5, a2=0.2, f1=80, Q=2e-8, f=80, b=5, mu0=1e-5, theta=0.08),
                         formulas=list(at="a1*t+a2", f1t="f1", Qt="Q*exp(theta*t)", ft="f", bt="b", mu0t="mu0*exp(theta*t)"), verbose=T) {
  # Optimization:
  res = optimize(data, start, formulas, verbose)
  invisible(res)
}