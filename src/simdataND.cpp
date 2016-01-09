// Simulation function for discrete case

/*******************************R-callable function***************************************/
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h> 
#include <iostream>
#include <RcppArmadillo.h>
//#include <random> - only for C++ 11
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

double mu(double t, double mu0, arma::mat b, arma::mat Q, double theta, arma::mat y);

double mu(double t, double mu0, arma::mat b, arma::mat Q, double theta, arma::mat y) {
    arma::mat res;
    res = (mu0 + y.t()*b + y.t()*Q*y)*exp(theta*t);
    //res = (mu0 + y.t()*b + y.t()*Q*y)*exp(theta);
    return res(0,0);
}

RcppExport SEXP simdata_ND(SEXP n, SEXP u_, SEXP R_, SEXP epsilon_, SEXP mu0_, SEXP b_, SEXP Q_, SEXP theta_, SEXP tstart_, SEXP ystart_, SEXP tend_, SEXP k_, SEXP dt_) {
    long N = as<long>(n); // Number of individuals
    int k = as<long>(k_); // Number of dimensions
    arma::mat u = as<arma::mat>(u_);
    //cout << u << endl;
    arma::mat R = as<arma::mat>(R_);
    //cout << R << endl;
    arma::mat epsilon = as<arma::mat>(epsilon_);
    //cout << epsilon << endl;
    double mu0 = as<double>(mu0_);
    arma::mat b = as<arma::mat>(b_);
    //cout << b << endl;
    arma::mat Q = as<arma::mat>(Q_);
    //cout << Q << endl;
    double theta  = as<double>(theta_);
    double tstart  = as<double>(tstart_);
    arma::mat ystart = as<arma::mat>(ystart_);
    
    double tend  = as<double>(tend_);
    //End of data loading
    double id = 0;
    double t1;
    double t2;
    double dt=as<double>(dt_);
    arma::mat y1;
    arma::mat y2;
    bool new_person = false;
    double S; 
    
    /*
    std::random_device rd;
    std::default_random_engine generator(rd()); 
    //std::default_random_engine generator(rand() % 100);
    std::uniform_real_distribution<double> uni_distr(0.0,1.0);
    std::uniform_real_distribution<double> uni_distr_dt(0.0,dt);
    */
    
    std::vector< std::vector<double> > data;
        
    for(int i=0; i < N; i++) { // Across all individuals
      t1 = tstart;
      t2 = t1 + dt;
      y1 = ystart;
      
      new_person = false;
      id = id + 1;
      
      while(new_person == false) {
        S = exp(-1*dt*mu(t1, mu0, b, Q, theta, y1));
        
        double xi = 0; // case (0 - alive, 1 - dead) indicator
        
        if(S > Rcpp::runif(1, 0.0, 1.0)[0]) { //if(S > uni_distr(generator)) {
          xi = 0; // case (0 - alive, 1 - dead) indicator
          arma::mat eps(k,1);
          
          for(int ii=0; ii < k; ii++) {
            /*
            std::normal_distribution<double> norm_distr(0.0,epsilon(ii,0));
            eps(ii,0) = norm_distr(generator);
            */
            eps(ii,0) = Rcpp::rnorm(1, 0.0, epsilon(ii,0.0))[0];
          }
          
          y2 = u + R * y1 + eps;
          new_person = false;
          
        } 
        else {
          xi = 1;
          y2 = arma::mat(k,1);
          /*t2 = t1 + uni_distr_dt(generator);*/
          t2 = t1 + Rcpp::runif(1, 0.0, dt)[0];
          new_person = true;
        }
          
        std::vector<double> row; row.resize(4+2*k);
        row[0] = id; row[1] = xi; row[2] = t1; row[3] = t2;
        
        int j=0;
        for(int ii=4; ii<(4 + 2*k-1); ii+=2) {
            row[ii] = y1(j,0);
            row[ii+1] = y2(j,0);
            j += 1;
        }
        
        data.push_back(row);
          
        if (new_person == false) {
          t1 = t2;
          t2 = t1 + dt;
        
          if(t2 > tend+dt) {
            new_person = true;
            break;
          }
        
          y1 = y2;
        }
     }
  }
  
  arma::mat res(data.size(), 4+2*k);
  for(int i=0; i<data.size(); i++) {
    for(int j=0; j<4+2*k; j++) {
      res(i, j) = data[i][j];
    }
  }
  
  return(Rcpp::wrap(res));

}