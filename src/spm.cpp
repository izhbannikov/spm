/*******************************R-callable function***************************************/
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

/*============GLOBAL VARIABLES===========*/
/*=====END OF GLOBAL VARIABLES===========*/

/*============FUNCTION DEFINITIONS*/
double mu(double t, double y1, double gamma1, double fH, double f1H, double mu0H, double thetaH, double QH, double *par);
double* func1(double t, double *y, double fH, double f1H, double aH, double bH, double QH);
double Q(double t, double QH);
/*=========END OF FUNCTION DEFINITIONS*/

double mu(double t, double y1, double gamma1, double fH, double f1H, double mu0H, double thetaH, double QH) {
  double hfH, hf1H, mu0Ht,  mu;
    
  hfH = fH-y1;
  hf1H = f1H-y1;  
  
  mu0Ht = mu0H*exp(thetaH*t);
  mu = mu0Ht + pow(hfH,2)*QH + QH*gamma1;
  
  return mu;
}

//Calculating m (y[1]) & gamma(y[2]):
double* func1(double t, double *y, double fH, double f1H, double aH, double bH, double QH) {
  double hfH, hf1H, dy1, dy2;
  hfH = fH-y[1];
  hf1H = f1H-y[1];
  dy1 = -1.0*aH*hf1H + 2.0*y[2]*Q(t, QH)*hfH;
  dy2 = 2.0*aH*y[2] + bH - 2.0*pow(y[2],2)*Q(t, QH); //dy2 <- 2*aH*y[2] + bH^2 - 2*y[2]^2*Q(t);
  double *res = new double[2];
  res[0] = dy1;
  res[1] = dy2;
  
  return res;
}

double Q(double t, double QH) {
  double Q;
  Q = QH;
  return Q;
}

RcppExport SEXP complik(SEXP dat, SEXP n, SEXP m, SEXP ah, SEXP f1h, SEXP qh, SEXP bh, SEXP fh, SEXP mu0h, SEXP thetah) {
    
    long N = as<long>(n); //Number of rows
    long M = as<long>(m); // Number of columns
    double aH = as<double>(ah);
    double f1H = as<double>(f1h);
    double QH = as<double>(qh);
    double bH = as<double>(bh);
    double fH = as<double>(fh);
    double mu0H = as<double>(mu0h); 
    double thetaH  = as<double>(thetah);
    //Actual data set
    Rcpp::NumericMatrix dd = Rcpp::NumericMatrix(dat);   
    
    /*End of data loading*/
    double *out, *k1ar, *yfin, *ytmp, *k2ar, *k3ar, *k4ar;
    out = new double[2];
    k1ar = new double[2];
    yfin = new double[2];
    ytmp = new double[2];
    k2ar = new double[2];
    k3ar = new double[2];
    k4ar = new double[2];
      
    double L; // Likelihood
    for(int i=0; i<N; i++) {
      //Solving differential equations on intervals:
      double t1 = dd(i,1); 
      double t2 = dd(i,2);
      double y1 = dd(i,3);
      double y2 = dd(i,4);
  
      int  nsteps = 2;
      double h=(t2-t1)/nsteps;
    
      //Integration:
      double s = h/3*(-1)*mu(t1,y1,0, fH, f1H, mu0H, thetaH, QH);
      double t = t1;
      out[0] = y1;
      out[1] = 0;
      int ifactor;
      
      for(int j = 0; j < nsteps; j++) {
         //Runge-Kutta method:
         k1ar = func1(t,out, fH, f1H, aH, bH, QH);
         yfin[0] = out[0] + h/6.00*k1ar[0];
         yfin[1] = out[1] + h/6.00*k1ar[1];
         ytmp[0] = out[0] + h/2.00*k1ar[0];
         ytmp[1] = out[1] + h/2.00*k1ar[1];
         
         k2ar = func1(t,ytmp, fH, f1H, aH, bH, QH);
         yfin[0] = yfin[0] + h/3.00*k2ar[0];
         yfin[1] = yfin[1] + h/3.00*k2ar[1];
         ytmp[0] = out[0] + h/2.00*k2ar[0];
         ytmp[1] = out[1] + h/2.00*k2ar[1];
        
         k3ar = func1(t,ytmp, fH, f1H, aH, bH, QH);
         yfin[0] = yfin[0] + h/3.00*k3ar[0];
         yfin[1] = yfin[1] + h/3.00*k3ar[1];
         ytmp[0] = out[0] + h*k3ar[0];
         ytmp[1] = out[1] + h*k3ar[1];
        
         k4ar = func1(t,ytmp, fH, f1H, aH, bH, QH);
         out[0] = yfin[0] + h/6.00*k4ar[0];
         out[1] = yfin[1] + h/6.00*k4ar[1];
         
         t = t + h;
      
        //Integration:
        if (j == nsteps) {
          ifactor = 1;
        } else {
          if ((j % 2) == 0) {
            ifactor = 2;
          } else {
            ifactor = 4;
          }
        }
        s = s + ifactor*h/3.00*(-1)*mu(t,y1,0, fH, f1H, mu0H, thetaH, QH);
      
      }
      //std::cout << s << "\n";
      
      double m2 = out[1];
      double gamma2 = out[2];
      double pi = 3.141592654;
    
      if(dd(i,0) == 0) { 
        double exp = -0.5*log(2*pi*gamma2)-pow((m2-y2),2)/2/gamma2;
        L = L + s + exp;
      } else {
        double logprobi = log(1 - exp(-1*mu(t2, m2, gamma2, fH, f1H, mu0H, thetaH, QH)));
        L = L + s + logprobi;
      }
    }
    //std::cout << L << "\n";
    return(Rcpp::wrap(L));
}
