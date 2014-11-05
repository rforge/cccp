/*
 * Inner product of two vectors in S.
 * sdot_nls is used for nonlinear, linear and second-order cone constraints.
 * sdot_p is used for positive semi-definite constraints.
*/

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(".sdot_nls")]]
double sdot_nls(SEXP us, SEXP vs) {

  Rcpp::NumericVector u(us);
  Rcpp::NumericVector v(vs);
  arma::colvec ua(u.begin(), u.size(), false);
  arma::colvec va(v.begin(), v.size(), false);

  double ans = arma::dot(ua, va);

  return ans;
}

// [[Rcpp::export(".sdot_p")]]
double sdot_p(SEXP us, SEXP vs, SEXP ms) {

  Rcpp::NumericVector u(us);
  Rcpp::NumericVector v(vs);
  int m = Rcpp::as<int>(ms);
  int n = u.size();
  double ans = 0.0;

  // squaring and summing diagonal elements
  for(int i = 0; i < n; i += m + 1){
    ans += u[i] * v[i];
  }
  // product of lower-diagonal elements, multiplied by two
  for(int i = 0; i < m; i++){
    for(int j = 0; j < (m - 1); j++){
      if(j < i){
	ans += 2.0 * u[i + m * j] * v[i + m * j];
      }
    }
  }

  return ans;
}
