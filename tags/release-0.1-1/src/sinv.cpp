/*
 * The inverse product x := (y 0\ x)
 */

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include "cccp.h"

// [[Rcpp::export(".sinv_nl")]]
SEXP sinv_nl(SEXP us, SEXP vs) {

  Rcpp::NumericMatrix u(us);
  Rcpp::NumericMatrix v(vs);
  Rcpp::NumericMatrix a(Rcpp::clone(us));

  int n = u.nrow();
  for(int i = 0; i < n; i++){
    a(i, 0) = u(i, 0) / v(i, 0);
  }

  return(Rcpp::wrap(a));
}

// [[Rcpp::export(".sinv_s")]]
SEXP sinv_s(SEXP us, SEXP vs) {

  Rcpp::NumericMatrix u(us);
  Rcpp::NumericMatrix v(vs);
  Rcpp::NumericMatrix a(Rcpp::clone(us));
  int n = u.nrow();

  double aa = jdot(v, v);
  double cc = u(0, 0);
  double dd = 0.0;

  for(int i = 1; i < n; i++){
    dd += u(i, 0) * v(i, 0);
  }
  a(0, 0) = cc * v(0, 0) - dd;
  for(int i = 1; i < n; i++){
    a(i, 0) = aa / v(0, 0) * u(i, 0);
    a(i, 0) += (dd / v(0, 0) - cc) * v(i, 0);
  }
  for(int i = 0; i < n; i++){
    a(i, 0) /= aa;
  }

  return(Rcpp::wrap(a));
}

// [[Rcpp::export(".sinv_p")]]
SEXP sinv_p(SEXP us, SEXP vs, SEXP ms) {

  arma::mat u = Rcpp::as<arma::mat>(us);
  arma::mat v = Rcpp::as<arma::mat>(vs);
  int m = Rcpp::as<int>(ms);
  u.reshape(m, m);
  v.reshape(m, m);
  arma::mat a = u;

  for(int i = 0; i < m; i++){
    for(int j = 0; j < m; j++){
      a.at(i,j) = u.at(i,j) * 2.0 / (v.at(i,i) + v.at(j,j));
	}
  }
  a.reshape(m * m, 1);

  return(Rcpp::wrap(a));
}
