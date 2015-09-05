#'Stochastic Process Modelling (SPM)
#'

spm <- function(dat,k=2, verbose=F, tol=NULL) {
  # Main function for Stochastic Process Modelling package
  # Parameters: 
  # - dat - input data frame
  # - k - number of dimensions
  # - tol - tolerance threshold for matrix inversion
  #
  # Step 1: estimation of starting point with quick discrete optimization:
  pars=spm_quick_MD(dat=dat[[2]],k=k)
  cat("Starting parameters:\n")
  
  # Step2: continuous slow optimization:
  data <- dat[[1]][,2:dim(dat[[1]])[2]]
  
  if(det(pars$pars2$Q) < 0) {
    cat("Error: determinant of Q < 0\n")
    cat("Q:\n")
    print(pars$pars2$Q)
    cat("Det(Q):\n")
    print(det(pars$pars2$Q))
    stop()
  }
  
  
  dd <- new.env()
  spm_integral_MD(data, c(pars$pars2$a, pars$pars2$f1, pars$pars2$Q, pars$pars2$b, pars$pars2$f, pars$pars2$mu0, pars$pars2$theta), k, dd, verbose)
  print(dd$results)
  res=list(starting=list(QH=pars$pars2$Q, aH=pars$pars2$a, bH=pars$pars2$b, f1H=pars$pars2$f1, fH=pars$pars2$f, mu0H=pars$pars2$mu0, theta=pars$pars2$theta), 
           estimated=dd$results)
  invisible(res)
}