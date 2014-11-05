#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export(".ssnt_l")]]
SEXP ssnt_l(SEXP us, SEXP Ws, SEXP invs) {

  Rcpp::NumericMatrix u(Rcpp::clone(us));
  Rcpp::List W(Ws);
  bool inv = Rcpp::as<bool>(invs);
  Rcpp::NumericMatrix w = Rcpp::as<Rcpp::NumericMatrix>(W["d"]);
  if(inv){
    w = Rcpp::as<Rcpp::NumericMatrix>(W["di"]);
  }

  int m = u.nrow();
  int n = u.ncol();
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      u(i, j) *= w(i, 0);
    }
  }

  return(Rcpp::wrap(u));
}

// [[Rcpp::export(".ssnt_s")]]
SEXP ssnt_s(SEXP us, SEXP Ws, SEXP invs) {

  arma::mat u = Rcpp::as<arma::mat>(Rcpp::clone(us));
  Rcpp::List W(Ws);
  arma::mat v = Rcpp::as<arma::mat>(W["v"]);
  double beta = Rcpp::as<double>(W["beta"]);
  bool inv = Rcpp::as<bool>(invs);
  double a = 0.0;
  arma::mat w;

  if(inv){
    u.row(0) = -u.row(0);
  }
  w = u.t() * v;
  u.row(0) = -u.row(0);
  u = 2 * v * w.t() + u;
  if(inv){
    u.row(0) = -u.row(0);
    a = 1 / beta;
  } else {
    a = beta;
  }
  u = a * u;

  return(Rcpp::wrap(u));
}


