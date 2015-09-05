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
  print(pars)
  # Step2: continuous slow optimization:
  data <- dat[[1]][,2:dim(dat[[1]])[2]]
  QH = pars$Q
  if(det(QH) < 0) {
    cat("Error: determinant of Q < 0\n")
    cat("Q:\n")
    print(QH)
    cat("Det(Q):\n")
    print(det(QH))
    stop()
  }
  aH = pars$R - 1
  bH = as.matrix(pars$eps, nrow=1)
  f1H = (-1)*pars$u %*% solve(aH, tol=tol)
  fH = -0.5 * pars$b %*% solve(QH, tol=tol)
  mu0H = pars$mu0 - fH %*% QH %*% t(fH)
  thetaH = pars$theta
  print(c(aH, f1H, QH, bH, fH, mu0H, thetaH,k))
  dd <- new.env()
  ans <- spm_integral_MD(data, c(aH, f1H, QH, bH, fH, mu0H, thetaH), k, dd)
  print(ans)
  res=list(starting=list(QH=QH, aH=aH, bH=bH, f1H=f1H, fH=fH, mu0H=mu0H, thetaH=thetaH), 
           estimated=ans)
  invisible(res)
}