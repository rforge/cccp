#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(".sdot_nls")]]
double sdot_nls(const arma::colvec & u, const arma::colvec & v) {
  double ans = arma::dot(u, v);
  return ans;
}
