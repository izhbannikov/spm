library(deSolve)

Q <- function(t) {
  Q <- QH*exp(thetaH*t)
  Q
}

func1 <- function(t, y, parms) {
  # y[[1]] = m, y[[2]] = gamma
  aH=parms$aH
  f1H=parms$f1H
  fH=parms$fH
  bH=parms$bH
  y1 = matrix(y[1:2],nrow=2,ncol=1, byrow=T)
  y2 = matrix(y[3:6], nrow=2, ncol=2, byrow=T)
  hfH <- t(fH) - y1
  hf1H <- t(f1H) - y1
  
  dm <- -1.0*aH%*%hf1H + 2.0*y2%*%Q(t)%*%hfH
  dgamma <- aH%*%y2 + y2%*%t(aH) + bH%*%t(bH) - 2.0*y2%*%Q(t)%*%y2 
  
  return(list(c(dm=dm, dgamma=dgamma)))
}

y0 = matrix(c(82.62224167,30.64715285), nrow=2, ncol=1, byrow=T)
gamma0 = matrix(0, nrow=2, ncol=2, byrow=T)

parms = list(aH=aH, f1H=f1H, QH=QH, fH=fH, bH=bH)

out <- ode(y = c(dm=y0, dgamma=gamma0), times = seq(50.19746171, 52.12581782, by=0.1),
           func = func1, parms = parms)

