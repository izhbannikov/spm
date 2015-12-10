#'Stochastic Process Modelling (SPM)
#'A main function that estimates parameters a, f1, Q, f, b, mu0, theta
#'from given dataset.
#'@param dat A dataset.
#'@param k Number of dimensions.
#'@param verbose A verbosing output indicator.
#'@param tol A tolerance threshold for matrix inversion.
#'@return A list of (1) Estimated starting point (from quick discrete optimization) and 
#'(2) Estimated coefficients.
#' @examples
#'library(spm)
#'#Prepare data for optimization
#'longdat <- read.csv(system.file("data","longdat.csv",package="spm"))
#'vitstat <- read.csv(system.file("data","vitstat.csv",package="spm"))
#' data=prepare_data(longdat=longdat, vitstat=vitstat,interval=1, col.status="IsDead", col.id="ID", col.age="Age", col.age.event="LSmort", covariates=c("DBP"), verbose=T)
#'#Parameters estimation:
#'pars=spm(data,k = 1)
#'pars
spm <- function(dat,k=2, verbose=F, tol=NULL) {
  # Main function for Stochastic Process Modelling package
  # Parameters: 
  # - dat - input data frame
  # - k - number of dimensions
  # - tol - tolerance threshold for matrix inversion
  #
  # Step 1: estimation of starting point with quick discrete optimization:
  pars=spm_discrete_MD(dat=dat[[2]],k=k)
  
  if(verbose) {
    cat("Starting parameters:\n")
    print(pars)
  }
  
  # Step2: continuous slow optimization:
  data <- dat[[1]][,2:dim(dat[[1]])[2]]
  
  if(det(pars$pars2$Q) < 0) {
    cat("Error: determinant of Q < 0\n")
    cat("Q:\n")
    print(pars$pars2$Q)
    cat("Det(Q):\n")
    print(det(pars$pars2$Q))
    res=list(starting=list(Q=pars$pars2$Q, a=pars$pars2$a, b=pars$pars2$b, f1=pars$pars2$f1, f=pars$pars2$f, mu0=pars$pars2$mu0, theta=pars$pars2$theta), 
             estimated=NA)
    
  } else {
  
    spm_continuous_MD(data, 
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