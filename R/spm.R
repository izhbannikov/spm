#'Stochastic Process Modelling (SPM)
#'A main function that estimates parameters a, f1, Q, f, b, mu0, theta
#'from given dataset.
#'@param x A dataset.
#'@param model A model type: c("continuous", "time-dependent").
#'@param formulas A list of formulas used in time-dependent model.
#'@param verbose A verbosing output indicator.
#'@param tol A tolerance threshold for matrix inversion.
#'@return A list of (1) coefficient estimates from discrete model and 
#'(if appropriate model type was provided): 
#'(2) coefficient estimates from continuous model and
#'(3) coefficient estimates from time-dependent model.
#' @examples
#'library(spm)
#'#Prepare data for optimization
#'longdat <- read.csv(system.file("data","longdat.csv",package="spm"))
#'vitstat <- read.csv(system.file("data","vitstat.csv",package="spm"))
#' data=prepare_data(longdat=longdat, vitstat=vitstat,interval=1, col.status="IsDead", col.id="ID", col.age="Age", col.age.event="LSmort", covariates=c("DBP"), verbose=T)
#'#Parameters estimation:
#'pars=spm(data)
#'pars
spm <- function(x, model=NULL, formulas = NULL, verbose=F, tol=NULL) {
  # Main function for Stochastic Process Modelling package
  # Parameters: 
  # - dat - input data frame
  # - tol - tolerance threshold for matrix inversion
  
  # List of available models:
  models <- c("continuous", "time-dependent")
  
  if(!is.null(model)) {
    if(length(unique(model)) > 2) {
      stop("Only continuous and/or time-dependent models can be used.")
    }
    for(i in 1:length(model)) {
      if(!(model[i] %in% models)) {
        cat("Warning: ", model[i], " - unknown model type!")
      }
    }
  }
  
  # Number of variables (dimensions):
  k <- (dim(x[[1]])[2] - 4)/2
  
  # Estimation of starting point with discrete optimization:
  pars=spm_discrete(dat=x[[2]],k=k)
  
  if(verbose) {
    cat("Starting parameters:\n")
    print(pars)
  }
  
  res <- list(estimated.discrete=list(Q=pars$pars2$Q, a=pars$pars2$a, b=pars$pars2$b, f1=pars$pars2$f1, f=pars$pars2$f, mu0=pars$pars2$mu0, theta=pars$pars2$theta))
  
  #If required, continuous-time optimization:
  if("continuous" %in% model) {
    data <- x[[1]][,2:dim(x[[1]])[2]]
  
    if(det(pars$pars2$Q) < 0) {
      cat("Error: determinant of Q < 0\n")
      cat("Q:\n")
      print(pars$pars2$Q)
      cat("Det(Q):\n")
      print(det(pars$pars2$Q))
      
      res[["estimated.continuous"]] <- NA
    
    } else {
  
      spm_continuous(data, 
                    a=pars$pars2$a, 
                    f1=pars$pars2$f1, 
                    Q=pars$pars2$Q, 
                    f=pars$pars2$f, 
                    b=pars$pars2$b, 
                    mu0=pars$pars2$mu0, 
                    theta=pars$pars2$theta, 
                    k, 
                    verbose)
 
      res[["estimated.continuous"]] <- get("results",envir=.GlobalEnv)
    }
  }
  
  if("time-dependent" %in% model) {
    data <- x[[1]][,2:dim(x[[1]])[2]]
    if(k > 1) {
      stop("Number of variables > 1. Model with time-dependent parameters can be used only with one variable!")
    }
    
    spm_time_dep <- function(data, 
                             start=list(a1=0, a2=estimated.discrete$a, f1=80, Q=2e-8, f=80, b=5, mu0=1e-5, theta=0.08),
                             formulas=list(at="a1*t+a2", f1t="f1", Qt="Q*exp(theta*t)", ft="f", bt="b", mu0t="mu0*exp(theta*t)"), 
                             verbose=TRUE,
                             lower_bound=NULL, upper_bound=NULL, factr=1e-16, lmult=0.5, umult=2) {
      
    
    spm_continuous(data, 
                     a=pars$pars2$a, 
                     f1=pars$pars2$f1, 
                     Q=pars$pars2$Q, 
                     f=pars$pars2$f, 
                     b=pars$pars2$b, 
                     mu0=pars$pars2$mu0, 
                     theta=pars$pars2$theta, 
                     k, 
                     verbose)
      
      res=list(starting=list(Q=pars$pars2$Q, a=pars$pars2$a, b=pars$pars2$b, f1=pars$pars2$f1, f=pars$pars2$f, mu0=pars$pars2$mu0, theta=pars$pars2$theta), 
               estimated=get("results",envir=.GlobalEnv))
    
  }
  
  invisible(res)
}