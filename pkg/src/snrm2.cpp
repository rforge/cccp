#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

double sdot_nls(const arma::colvec & u, const arma::colvec & v);
double sdot_p(const arma::colvec & u, const arma::colvec & v, const int & m);

/*
Norm of a vector in S.
snrm2_nls is used for nonlinear, linear and second-order cone constraints.
snrm2_p is used for positive semi-definite constraints.
*/

// [[Rcpp::export(".snrm2_nls")]]
double snrm2_nls(const arma::colvec & u) {
  double ans = sqrt(sdot_nls(u, u));
  return ans;
}

// [[Rcpp::export(".snrm2_p")]]
double snrm2_p(const arma::colvec & u, const int & m) {
  double ans = 0.0;
  ans = sqrt(sdot_p(u, u, m));
  return ans;
}
