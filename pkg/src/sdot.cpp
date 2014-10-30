#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

/*
Inner product of two vectors in S.
sdot_nls is used for nonlinear, linear and second-order cone constraints.
sdot_p is used for positive semi-definite constraints.
*/


// [[Rcpp::export(".sdot_nls")]]
double sdot_nls(const arma::colvec & u, const arma::colvec & v) {
  double ans = arma::dot(u, v);
  return ans;
}

// [[Rcpp::export(".sdot_p")]]
double sdot_p(const arma::colvec & u, const arma::colvec & v, const int & m) {
  double ans = 0.0;
  int n = u.n_elem;
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
