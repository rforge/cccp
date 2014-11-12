/*
 * Returns x' * J * y, whereby J = [1, 0; 0, -I]
*/

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(".jdot")]]
double jdot(SEXP us, SEXP vs) {

  arma::mat ua = Rcpp::as<arma::mat>(us);
  arma::mat va = Rcpp::as<arma::mat>(vs);
  double a = 0.0;
  int n = ua.n_rows;

  a = ua(0, 0) * va(0, 0);
  for(int i = 1; i < n; i++){
    a -= ua(i, 0) * va(i, 0);
  }
  return(a);
}

// [[Rcpp::export(".jnrm2")]]
double jnrm2(SEXP us) {
  double a = sqrt(jdot(us, us));
  return(a);
}
