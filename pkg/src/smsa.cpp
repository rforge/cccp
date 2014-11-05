/*
 * Applying maximum step length
*/

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(".smsa_nl")]]
SEXP smsa_nl(SEXP us, SEXP alphas, SEXP inits) {
  Rcpp::NumericMatrix u(us);
  int n = u.nrow();
  double alpha = Rcpp::as<double>(alphas);
  bool init = Rcpp::as<bool>(inits);

  if(init){
    for(int i = 0; i < n; i++){
      u(i, 0) = u(i, 0) + (1 + alpha);
    }
  } else {
    for(int i = 0; i < n; i++){
      u(i, 0) = 1.0 + alpha * u(i, 0);
    }
  }

  return(Rcpp::wrap(u));
}

// [[Rcpp::export(".smsa_s")]]
SEXP smsa_s(SEXP us, SEXP alphas, SEXP inits) {
  Rcpp::NumericMatrix u(us);
  int n = u.nrow();
  double alpha = Rcpp::as<double>(alphas);
  bool init = Rcpp::as<bool>(inits);

  if(init){
    u(0, 0) += 1.0 + alpha;
  } else {
    for(int i = 0; i < n; i++){
      u(i, 0) = alpha * u(i, 0);
    }
    u(0, 0) += 1.0;
  }

  return(Rcpp::wrap(u));
}

// [[Rcpp::export(".smsa1_p")]]
SEXP smsa1_p(SEXP us, SEXP alphas, SEXP ms) {

  arma::mat u = Rcpp::as<arma::mat>(us);
  double alpha = Rcpp::as<double>(alphas);
  int m = Rcpp::as<int>(ms);

  u.reshape(m, m);
  u.diag() = u.diag() + (1 + alpha);
  u.reshape(m * m, 1);

  return(Rcpp::wrap(u));
}

// [[Rcpp::export(".smsa2_p")]]
SEXP smsa2_p(SEXP us, SEXP alphas, SEXP sigmas, SEXP lambdas, SEXP ms) {

  arma::mat u = Rcpp::as<arma::mat>(us);
  double alpha = Rcpp::as<double>(alphas);
  Rcpp::NumericVector s(Rcpp::clone(sigmas));
  arma::mat l = Rcpp::as<arma::mat>(lambdas);
  int m = Rcpp::as<int>(ms);

  u.reshape(m, m);
  l.reshape(m, m);
  for(int i = 0; i < m; i++){
    s(i) = 1 + alpha * s(i);
    u.col(i) = u.col(i) * sqrt(s(i) / l(i, i)); 
  }
  u.reshape(m * m, 1);

  return(Rcpp::wrap(u));
}
