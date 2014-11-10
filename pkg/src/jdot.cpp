/*
 * Returns x' * J * y, whereby J = [1, 0; 0, -I]
*/

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(".jdot")]]
double jdot(SEXP us, SEXP vs) {

  Rcpp::NumericMatrix u(us);
  Rcpp::NumericMatrix v(vs);
  double a = 0.0;
  int n = u.nrow();

  a = u(0, 0) * v(0, 0); 
  for(int i = 1; i < n; i++){
    a -= u(i, 0) * v(i, 0);
  }

  return(a);
}

// [[Rcpp::export(".jnrm2")]]
double jnrm2(SEXP us) {
  double a = sqrt(jdot(us, us));
  return(a);
}
