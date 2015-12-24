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
spm <- function(x, model="discrete", formulas = NULL, verbose=F, tol=NULL) {
  
  # List of available models:
  models <- c("discrete", "continuous", "time-dependent")
  
  if(!(model %in% models)) {
    stop(cat(model, " - unknown model type!"))
  }
  
  # Number of variables (dimensions):
  k <- (dim(x[[1]])[2] - 4)/2
  
  
  if(model == "discrete") {
    # Estimation of starting point with discrete optimization:
    pars <- spm_discrete(dat=x[[2]],k=k)
    res <- list(Ak2005=list(u=pars$pars1$u, 
                            R=pars$pars1$R, 
                            b=pars$pars1$b, 
                            Q=pars$pars1$Q, 
                            epsilon=pars$pars1$eps,
                            mu0=pars$pars1$mu0,
                            theta=pars$pars1$theta), 
                Ya2007=list(a=pars$pars2$a, 
                            f1=pars$pars2$f1,
                            Q=pars$pars2$Q,
                            f=pars$pars2$f, 
                            b=pars$pars2$b, 
                            mu0=pars$pars2$mu0, 
                            theta=pars$pars2$theta))
    
  }
  
  
  if(model == "continuous") {
    pars <- spm_discrete(dat=x[[2]],k=k)
    data <- x[[1]][,2:dim(x[[1]])[2]]
  
    if(verbose) {
      cat("Starting parameters:\n")
      print(pars)
    }
    
    if(det(pars$pars2$Q) < 0) {
      cat("Error: determinant of Q < 0\n")
      cat("Q:\n")
      print(pars$pars2$Q)
      cat("Det(Q):\n")
      print(det(pars$pars2$Q))
      
      res <- NA
    
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
  
      res.t <- get("results",envir=.GlobalEnv)
      
      Q.c <- res.t$Q
      R.c <- res.t$a + diag(k)
      eps.c <- as.matrix(res.t$b)
      u.c <- (-1)*(res.t$f1 %*% res.t$a)
      b.c <- -2*res.t$f %*% res.t$Q
      mu0.c <- res.t$mu0 + res.t$f %*% res.t$Q %*% t(res.t$f)
      theta.c <- res.t$theta
      
      res <- list(Ak2005=list(u=u.c, 
                              R=R.c, 
                              b=b.c, 
                              Q=Q.c, 
                              epsilon=eps.c,
                              mu0=mu0.c,
                              theta=theta.c), 
                  Ya2007=list(a=res.t$a, 
                              f1=res.t$f1,
                              Q=res.t$Q,
                              f=res.t$f, 
                              b=res.t$b, 
                              mu0=res.t$mu0, 
                              theta=res.t$theta))
      
    }
  }
  
  #if(model == "time-dependent") {
  #  data <- x[[1]][,2:dim(x[[1]])[2]]
  #  if(k > 1) {
  #    stop("Number of variables > 1. Model with time-dependent parameters can be used only with one variable!")
  #  }
  #  
  #  spm_time_dep <- function(data, 
  #                           start=list(a1=0, a2=estimated.discrete$a, f1=80, Q=2e-8, f=80, b=5, mu0=1e-5, theta=0.08),
  #                           formulas=list(at="a1*t+a2", f1t="f1", Qt="Q*exp(theta*t)", ft="f", bt="b", mu0t="mu0*exp(theta*t)"), 
  #                           verbose=TRUE,
  #                           lower_bound=NULL, upper_bound=NULL, factr=1e-16, lmult=0.5, umult=2) {
  #    
  #  
  #  spm_continuous(data, 
  #                   a=pars$pars2$a, 
  #                   f1=pars$pars2$f1, 
  #                   Q=pars$pars2$Q, 
  #                   f=pars$pars2$f, 
  #                   b=pars$pars2$b, 
  #                   mu0=pars$pars2$mu0, 
  #                   theta=pars$pars2$theta, 
  #                   k, 
  #                   verbose)
  #    
  #    res <- list(starting=list(Q=pars$pars2$Q, a=pars$pars2$a, b=pars$pars2$b, f1=pars$pars2$f1, f=pars$pars2$f, mu0=pars$pars2$mu0, theta=pars$pars2$theta), 
  #             estimated=get("results",envir=.GlobalEnv))
  #  
  #}
  
  invisible(res)
}