#'Stochastic Process Modelling (SPM)
#'

spm <- function(dat,k=2, tol=NULL) {
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
  
  if(det(pars$pars2$QH) < 0) {
    cat("Error: determinant of Q < 0\n")
    cat("Q:\n")
    print(pars$pars2$QH)
    cat("Det(Q):\n")
    print(det(pars$pars2$QH))
    stop()
  }
  
  #aH = pars$R - 1
  #bH = as.matrix(pars$eps, nrow=1)
  #f1H = (-1)*pars$u %*% solve(aH, tol=tol)
  #fH = -0.5 * pars$b %*% solve(QH, tol=tol)
  #mu0H = pars$mu0 - fH %*% QH %*% t(fH)
  #thetaH = pars$theta
  
  #print(c(pars$pars2$aH, pars$pars2$f1H, pars$pars2$QH, pars$pars2$bH, pars$pars2$fH, pars$pars2$mu0H, pars$pars2$thetaH,k))
  
  dd <- new.env()
  ans <- spm_integral_MD(data, c(pars$pars2$aH, pars$pars2$f1H, pars$pars2$QH, pars$pars2$bH, pars$pars2$fH, pars$pars2$mu0H, pars$pars2$thetaH), k, dd)
  print(ans)
  res=list(starting=list(QH=pars$pars2$QH, aH=pars$pars2$aH, bH=pars$pars2$bH, f1H=pars$pars2$f1H, fH=pars$pars2$fH, mu0H=pars$pars2$mu0H, thetaH=pars$pars2$thetaH), 
           estimated=ans)
  invisible(res)
}