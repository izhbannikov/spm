
/*******************************R-callable function***************************************/
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>  
//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

/*============GLOBAL VARIABLES===========*/
/*=====END OF GLOBAL VARIABLES===========*/

/*============FUNCTION DEFINITIONS*/

double mu(double t, arma::mat y1, arma::mat gamma1, arma::mat fH, arma::mat f1H, double mu0H, double thetaH, arma::mat QH);
//double* func1(double t, double *y, double fH, double f1H, double aH, double bH, double QH, double theta);
arma::mat Q(double t, arma::mat QH, double theta);

// Matrix operations:
Rcpp::NumericMatrix add(Rcpp::NumericMatrix a, Rcpp::NumericMatrix b);
Rcpp::NumericMatrix muliplyByConstant(double constant, Rcpp::NumericMatrix a);
Rcpp::NumericMatrix multiplyMatrices(Rcpp::NumericMatrix a, Rcpp::NumericMatrix b);
Rcpp::NumericMatrix powerMatrix(Rcpp::NumericMatrix a, int p);

//=========END OF FUNCTION DEFINITIONS

Rcpp::NumericMatrix add(Rcpp::NumericMatrix a, Rcpp::NumericMatrix b) {
  if(a.ncol() != b.ncol() || a.nrow() != b.nrow()) {
    std::cout << "Error: not conformable arrays.\n";
    return R_NaReal;
  } 
  
  Rcpp::NumericMatrix c(a.nrow(),a.ncol());
  for(int i=0; i<a.nrow(); i++) {
    for(int j=0; j<a.ncol(); j++) {
      c(i,j) = a(i,j) + b(i,j);
    }
  }
  return c;
}

Rcpp::NumericMatrix muliplyByConstant(double constant, Rcpp::NumericMatrix a) {
  Rcpp::NumericMatrix c(a.nrow(),a.ncol());
  for(int i=0; i<a.nrow(); i++) {
    for(int j=0; j<a.ncol(); j++) {
      c(i,j) = constant*a(i,j);
    }
  }
  return c;
}

Rcpp::NumericMatrix multiplyMatrices(Rcpp::NumericMatrix a, Rcpp::NumericMatrix b) {
  /*Function to multiply matrices*/
  if(a.ncol() != b.nrow()) {
    std::cout << "Error: not conformable arrays.\n";
    return R_NaReal;
  }
  
  Rcpp::NumericMatrix c(a.nrow(),b.ncol());
  for(int i=0; i<a.nrow(); i++) {
    for(int j=0; j<b.ncol(); j++) {
      for(int k=0; k<a.ncol(); k++) {
        c(i,j) += a(i,k) * b(k,j);
      }
    }
  }
  return c;
}

Rcpp::NumericMatrix powerMatrix(Rcpp::NumericMatrix a, int p) {
  Rcpp::NumericMatrix c(a.nrow(),a.ncol());
  for(int i=0; i<a.nrow(); i++) {
    for(int j=0; j<a.ncol(); j++) {
      c(i,j) = pow(a(i,j),p);
    }
  }
  return c;
}



double mu(double t, arma::mat y1, arma::mat gamma1, arma::mat fH, arma::mat f1H, double mu0H, double thetaH, arma::mat QH) {
  arma::mat hfH;
  arma::mat hf1H;
  double mu0Ht;
  arma::mat mu;
    
  hfH = fH - y1; //fH-y1;  
  hf1H = f1H - y1; //f1H-y1;  
  
  mu0Ht = mu0H*exp(thetaH*t);
  mu = mu0Ht + hfH.t()*QH*hfH + QH*gamma1;
  
  return mu(1,1);
}

/*
//Calculating m (y[1]) & gamma(y[2]):
double func1(double t, arma::mat y, arma::mat gamma, arma::mat fH, arma::mat f1H, arma::mat aH, arma::mat bH, arma::mat QH, double theta) {
  arma::mat hfH, hf1H, dy1, dy2;
  hfH = fH-y;
  hf1H = f1H-y;
  dy1 = -1.00*aH*hf1H + 2.00*gamma*Q(t, QH, theta)*hfH;
  dy2 = aH*gamma + gamma*aH.t() + bH*bH.t() - 2.00*gamma.t()*Q(t, QH, theta)*gamma; //dy2 <- 2*aH*y[2] + bH^2 - 2*y[2]^2*Q(t);
  double *res = new double[2];
  res[0] = dy1;
  res[1] = dy2;
  
  return res;
}
*/

arma::mat Q(double t, arma::mat QH, double theta) {
  arma::mat Q;
  Q = QH*exp(theta*t);
  return Q;
}


RcppExport SEXP complikMD(SEXP dat, SEXP n, SEXP m, SEXP ah, SEXP f1h, SEXP qh, SEXP bh, SEXP fh, SEXP mu0h, SEXP thetah, SEXP k) {
    
    long N = as<long>(n); // Number of rows
    long M = as<long>(m); // Number of columns
    arma::mat aH = as<arma::mat>(ah);
    arma::mat f1H = as<arma::mat>(f1h);
    arma::mat QH = as<arma::mat>(qh);
    arma::mat bH = as<arma::mat>(bh);
    arma::mat fH = as<arma::mat>(fh);
    double mu0H = as<double>(mu0h); 
    double thetaH  = as<double>(thetah);
    int dim = as<int>(k);
    //Actual data set
    arma::mat dd = as<arma::mat>(dat);  
    double L; // Likelihood
    
    //End of data loading
    arma::mat out, k1ar, yfin, ytmp, k2ar, k3ar, k4ar;
    /*
    L = 0;
    for(int i=0; i<N; i++) {
      //Solving differential equations on intervals:
      double t1 = dd(i,1); 
      double t2 = dd(i,2);
      //double y1 = dd(i,3);
      //double y2 = dd(i,4);
      arma::mat y1(dim,1); 
      arma::mat y2(dim,1);
      Rcpp::NumericMatrix gamma1(dim,1);
      for(int ii=0; ii<M; ii+=2) {
        y1(ii,1) = dd(i,ii+3);
        y2(ii,1) = dd(i,ii+3);
        gamma1(ii,1) = 0;
      }
      
      double  nsteps = 2.00;
      double h=(t2-t1)/nsteps;
      
      //Integration:
      double s = h/3.00*(-1.00)*mu(t1, y1, gamma1, fH, f1H, mu0H, thetaH, QH);
      double t = t1;
      out[0] = y1;
      out[1] = 0.00;
      double ifactor;
      
      for(int j = 0; j < nsteps; j++) {
         //Runge-Kutta method:
         k1ar = func1(t,out, fH, f1H, aH, bH, QH, thetaH);
         yfin[0] = out[0] + h/6.00*k1ar[0];
         yfin[1] = out[1] + h/6.00*k1ar[1];
         ytmp[0] = out[0] + h/2.00*k1ar[0];
         ytmp[1] = out[1] + h/2.00*k1ar[1];
         
         k2ar = func1(t,ytmp, fH, f1H, aH, bH, QH, thetaH);
         yfin[0] = yfin[0] + h/3.00*k2ar[0];
         yfin[1] = yfin[1] + h/3.00*k2ar[1];
         ytmp[0] = out[0] + h/2.00*k2ar[0];
         ytmp[1] = out[1] + h/2.00*k2ar[1];
        
         k3ar = func1(t,ytmp, fH, f1H, aH, bH, QH, thetaH);
         yfin[0] = yfin[0] + h/3.00*k3ar[0];
         yfin[1] = yfin[1] + h/3.00*k3ar[1];
         ytmp[0] = out[0] + h*k3ar[0];
         ytmp[1] = out[1] + h*k3ar[1];
        
         k4ar = func1(t,ytmp, fH, f1H, aH, bH, QH, thetaH);
         out[0] = yfin[0] + h/6.00*k4ar[0];
         out[1] = yfin[1] + h/6.00*k4ar[1];
         
         t = t + h;
      
        //Integration:
        if (j == nsteps-1) {
          ifactor = 1.00;
        } else {
          if (((j % 2) == 0) && (j != 0)) {
            ifactor = 2.00;
          } else {
            ifactor = 4.00;
          }
        }
        //cout << ifactor << "\n";
        s = s + ifactor*h/3.00*(-1.00)*mu(t,out[0],out[1], fH, f1H, mu0H, thetaH, QH);
        
      }
      
      double m2 = out[0];
      double gamma2 = out[1];
      double pi = 3.141592654;
    
      if(dd(i,0) == 0) { 
        double exp = -0.50*log(2.00*pi*gamma2)-pow((m2-y2),2.00)/2.00/gamma2;
        L = L + s + exp;
        //cout << exp << "\n";
      } else {
        double logprobi = log(1.00 - exp(-1.00*mu(t2, m2, gamma2, fH, f1H, mu0H, thetaH, QH)));
        L = L + s + logprobi;
        //cout << s << " " << logprobi << " " << m2 << " " << gamma2 << "\n";
      }
      //break;
    }
    //std::cout << L << "\n";
    return(Rcpp::wrap(L));
    */
    return(Rcpp::wrap(""));
}

RcppExport SEXP testMatrixMultiply(SEXP aa, SEXP bb) {
  Rcpp::NumericMatrix a = Rcpp::NumericMatrix(aa);
  Rcpp::NumericMatrix b = Rcpp::NumericMatrix(bb);
  
  return(Rcpp::wrap(multiplyMatrices(a, b)));
}

RcppExport SEXP testMatrixAdd(SEXP aa, SEXP bb) {
  Rcpp::NumericMatrix a = Rcpp::NumericMatrix(aa);
  Rcpp::NumericMatrix b = Rcpp::NumericMatrix(bb);
  
  return(Rcpp::wrap(add(a, b)));
}

RcppExport SEXP testMatrixMultiplyByConst(SEXP aa, SEXP bb) {
  double constant = as<double>(aa);
  Rcpp::NumericMatrix a = Rcpp::NumericMatrix(bb);
  
  return(Rcpp::wrap(muliplyByConstant(constant, a)));
}

RcppExport SEXP testPowerMatrix(SEXP aa, SEXP pp) {
  double p = as<double>(pp);
  Rcpp::NumericMatrix a = Rcpp::NumericMatrix(aa);
  
  return(Rcpp::wrap(powerMatrix(a,p)));
}




