
/*******************************R-callable function***************************************/
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h> 
#include <iostream>
//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

/*============GLOBAL VARIABLES===========*/
/*=====END OF GLOBAL VARIABLES===========*/

/*============FUNCTION DEFINITIONS*/

double mu(double t, arma::mat y1, arma::mat gamma1, arma::mat fH, arma::mat f1H, double mu0H, double thetaH, arma::mat QH);
void func1(arma::mat *res, double t, arma::mat *y, arma::mat fH, arma::mat f1H, arma::mat aH, arma::mat bH, arma::mat QH, double theta);
arma::mat Q(double t, arma::mat QH, double theta);


double mu(double t, arma::mat y1, arma::mat gamma1, arma::mat fH, arma::mat f1H, double mu0H, double thetaH, arma::mat QH) {
  arma::mat hfH;
  arma::mat hf1H;
  double mu0Ht;
  arma::mat mu;
  
  hfH = fH.t() - y1; //fH-y1;  
  hf1H = f1H.t() - y1; //f1H-y1;  
  mu0Ht = mu0H*exp(thetaH*t);
  arma::mat QH_gamma1 = QH*gamma1;
  mu = mu0Ht + (hfH.t()*QH)*hfH + arma::sum((QH_gamma1).diag());
  
  return mu(0,0);
}


//Calculating m (y[1]) & gamma(y[2]):
void func1(arma::mat *res, double t, arma::mat *y, arma::mat fH, arma::mat f1H, arma::mat aH, 
          arma::mat bH, arma::mat QH, double theta) {
  
  // y[0] = m, k by 1, y[1] = gamma, k by k, where k is # of dimensions
  // dy1 = dm/dt, dy2 = dgamma/dt
  
  arma::mat hfH, hf1H, dy1, dy2;
  hfH = fH.t() - y[0]; 
  hf1H = f1H.t() - y[0];
  dy1 = -1.00 * (aH*hf1H) + 2.00 * ((y[1]*Q(t, QH, theta))*hfH);
  dy2 = aH*y[1] + y[1]*aH.t() + bH*bH.t() - 2.00 * ((y[1]*Q(t, QH, theta))*y[1]);
  
  res[0] = dy1; res[1] = dy2;
  
}


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
    //std::cout << "1\n";
    arma::mat dd = as<arma::mat>(dat);  
    //std::cout << "passed\n";
    double L; // Likelihood
    
    //End of data loading
    arma::mat *out, *k1ar, *yfin, *ytmp, *k2ar, *k3ar, *k4ar;
    out = new arma::mat[2];
    k1ar = new arma::mat[2];
    yfin = new arma::mat[2];
    ytmp = new arma::mat[2];
    k2ar = new arma::mat[2];
    k3ar = new arma::mat[2];
    k4ar = new arma::mat[2];
    
    arma::mat y1(dim,1);
    arma::mat y2(dim,1);
    arma::mat gamma1(dim,dim);
    
    L = 0;
    double  nsteps = 2;
    //std::cout << N << std::endl;
    
    
    for(int i=0; i<N; i++) {
      //Solving differential equations on intervals:
      double t1 = dd(i,1); 
      double t2 = dd(i,2);
      //std::cout<< i << " " << t1 <<  " " << t2  << " " << M << std::endl;
    
      int jj=0;
      for(int ii=3; ii<M; ii+=2) {
        y1(jj,0) = dd(i,ii);
        y2(jj,0) = dd(i,ii+1);
        jj += 1;
      }
      
      double tdiff = t2-t1;
      if(tdiff > 2) {
        nsteps = 2*tdiff;
      }
      
      double h = tdiff/nsteps;
      
      //Integration:
      gamma1.zeros(); // set gamma1 to zero matrix
      //std::cout << "2\n";
      double s = h/3.00*(-1.00)*mu(t1, y1, gamma1, fH, f1H, mu0H, thetaH, QH);
      //std::cout << "passed 2\n";
      double t = t1;
      out[0] = y1;
      out[1] = gamma1;
      double ifactor;
      
      for(int j = 0; j < nsteps; j++) {
         //Runge-Kutta method:
         func1(k1ar, t, out, fH, f1H, aH, bH, QH, thetaH);
         yfin[0] = out[0] + h/6.00*k1ar[0];
         yfin[1] = out[1] + h/6.00*k1ar[1];
         ytmp[0] = out[0] + h/2.00*k1ar[0];
         ytmp[1] = out[1] + h/2.00*k1ar[1];
         if(isnan(yfin[1](0,0))) {
           cout << "k1ar:\n" << k1ar << "\nytmp:\n" << ytmp << "\n"; 
         }
         
         func1(k2ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
         yfin[0] = yfin[0] + h/3.00*k2ar[0];
         yfin[1] = yfin[1] + h/3.00*k2ar[1];
         ytmp[0] = out[0] + h/2.00*k2ar[0];
         ytmp[1] = out[1] + h/2.00*k2ar[1];
         if(isnan(yfin[1](0,0))) {
           cout << "k2ar:\n" << k2ar << "\nytmp:\n" << ytmp << "\n"; 
         }
         
         func1(k3ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
         yfin[0] = yfin[0] + h/3.00*k3ar[0];
         yfin[1] = yfin[1] + h/3.00*k3ar[1];
         ytmp[0] = out[0] + h*k3ar[0];
         ytmp[1] = out[1] + h*k3ar[1];
         if(isnan(yfin[1](0,0))) {
           cout << "k3ar:\n" << k3ar << "\nytmp:\n" << ytmp << "\n"; 
         }
         
         func1(k4ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
         out[0] = yfin[0] + h/6.00*k4ar[0];
         out[1] = yfin[1] + h/6.00*k4ar[1];
         if(isnan(out[1](0,0))) {
           cout << "k4ar:\n" << k4ar << "\nytmp:\n" << ytmp << "\n"; 
         }
         
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
      
      arma::mat m2 = out[0];
      arma::mat gamma2 = out[1];
      double pi = 3.141592654;
      
      if(dd(i,0) == 0) { 
        
        //cout << "QH:\n" << QH << "\n" << "gamma:\n" << gamma2 << "\nm:" << m2 << "\n" << y2 << "\nL:" << L << "\n";
        arma::mat exp = -0.50*dim*log(2.00*pi*det(gamma2)) - 0.50*(m2-y2).t()*pinv(gamma2,0.000000000001)*(m2-y2);
        //arma::mat exp = -0.50*dim*log(2.00*pi*det(gamma2)) - 0.50*(m2-y2).t()*inv(gamma2)*(m2-y2);
        L += s + exp(0,0);
        //cout << exp << endl;
        if((det(gamma2) < 0) && (det(QH) > 0)) {
          cout << "Det gamma < 0\n";
          cout << "i:\n" << i << "\nt1:\n" << t1 << "\nt2:\n" << t2 << "\ny1:\n" << y1 << "\ny2:\n" << y2 << " "<< m2 << "\n";
          cout << "QH:\n" << QH << ", det: " << det(QH) << "\ngamma:\n" << gamma2 << ", get: " << det(gamma2) << "\nbH:" << bH << "\n" << bH*bH.t() << "\nL:" << L << "\n";
          cout << nsteps << "\n";
          break;
        }
        if(isnan(exp(0,0))) {
          cout << "Det gamma < 0\n";
          cout << "i:\n" << i << "\nt1:\n" << t1 << "\nt2:\n" << t2 << "\ny1:\n" << y1 << "\ny2:\n" << y2 << "\n";
          cout << "QH:\n" << QH << ", det: " << det(QH) << "\ngamma:\n" << gamma2 << ", get: " << det(gamma2) << "\nbH:" << bH << "\n" << bH*bH.t() << "\nL:" << L << "\n";
          cout << nsteps << "\n";
          break;
        }
      } else {
        //cout << "??\n";
        double logprobi = log(1.00 - exp(-1.00*mu(t2, m2, gamma2, fH, f1H, mu0H, thetaH, QH)));
        L += s + logprobi;
        //cout << s << " " << logprobi << " " << m2 << " " << gamma2 << "\n";
      }
      
      
    }
    
    delete[] out;
    delete[] k1ar;
    delete[] yfin;
    delete[] ytmp;
    delete[] k2ar;
    delete[] k3ar;
    delete[] k4ar;
    
    //std::cout << L << "\n";
    return(Rcpp::wrap(L));
}

