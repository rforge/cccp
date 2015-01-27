/*
 * One element with respect to a cone
*/
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(".sone_nl")]]
SEXP sone_nl(SEXP us) {
  Rcpp::NumericMatrix u(Rcpp::clone(us));
  std::fill(u.begin(), u.end(), 1);
  return(Rcpp::wrap(u));
}

// [[Rcpp::export(".sone_s")]]
SEXP sone_s(SEXP us) {
  Rcpp::NumericMatrix u(Rcpp::clone(us));
  std::fill(u.begin(), u.end(), 0);
  u(0, 0) = 1.0;
  return(Rcpp::wrap(u));
}

// [[Rcpp::export(".sone_p")]]
SEXP sone_p(SEXP ms) {
  int m = Rcpp::as<int>(ms);
  arma::mat Ident = arma::eye(m, m);
  Ident.reshape(m * m, 1);
  return(Rcpp::wrap(Ident));
}
