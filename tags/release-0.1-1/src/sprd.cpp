/*
 * Returns x := (y o x)
*/

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(".sprd_nl")]]
SEXP sprd_nl(SEXP us, SEXP vs) {

  Rcpp::NumericMatrix u(us);
  Rcpp::NumericMatrix v(vs);
  Rcpp::NumericMatrix a(Rcpp::clone(us));
  int n = u.nrow();

  for(int i = 0; i < n; i++){
    a(i, 0) = u(i, 0) * v(i, 0);
  }

  return(Rcpp::wrap(a));
}

// [[Rcpp::export(".sprd_s")]]
SEXP sprd_s(SEXP us, SEXP vs) {

  arma::mat u = Rcpp::as<arma::mat>(us);
  arma::mat v = Rcpp::as<arma::mat>(vs);
  arma::mat a = v;
  int n = u.n_rows;

  a(0, 0) = arma::dot(u, v);
  for(int i = 1; i < n; i++){
    a(i, 0) = u(0, 0) * v(i, 0) + v(0, 0) * u(i, 0);
  }

  return(Rcpp::wrap(a));
}

// [[Rcpp::export(".sprd_p")]]
SEXP sprd_p(SEXP us, SEXP vs, SEXP ms) {

  int m = Rcpp::as<int>(ms);
  arma::mat u = Rcpp::as<arma::mat>(us);
  u.reshape(m, m);
  arma::mat v = Rcpp::as<arma::mat>(vs);
  v.reshape(m, m);
  arma::mat a(m, m);

  a = 0.5 * (u * v + v * u);
  a.reshape(m * m, 1);

  return(Rcpp::wrap(a));
}
