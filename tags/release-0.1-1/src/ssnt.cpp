#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(".ssnt_n")]]
SEXP ssnt_n(SEXP us, SEXP Ws, SEXP invs) {

  Rcpp::NumericMatrix u(Rcpp::clone(us));
  Rcpp::List W(Ws);
  bool inv = Rcpp::as<bool>(invs);
  arma::mat w;

  if(inv){
    w = Rcpp::as<arma::mat>(W["dnli"]);
  } else {
    w = Rcpp::as<arma::mat>(W["dnl"]);
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

// [[Rcpp::export(".ssnt_l")]]
SEXP ssnt_l(SEXP us, SEXP Ws, SEXP invs) {

  Rcpp::NumericMatrix u(Rcpp::clone(us));
  Rcpp::List W(Ws);
  bool inv = Rcpp::as<bool>(invs);
  arma::mat w;

  if(inv){
    w = Rcpp::as<arma::mat>(W["di"]);
  } else {
    w = Rcpp::as<arma::mat>(W["d"]);
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

// [[Rcpp::export(".ssnt_p")]]
SEXP ssnt_p(SEXP us, SEXP Ws, SEXP transs, SEXP invs) {

  arma::mat u = Rcpp::as<arma::mat>(Rcpp::clone(us));
  Rcpp::List W(Ws);
  bool inv = Rcpp::as<bool>(invs);
  bool trans = Rcpp::as<bool>(transs);
  bool tt;
  int n, m;
  arma::mat w, U, a, ans;

  if(inv){
    w = Rcpp::as<arma::mat>(W["rti"]);
    tt = trans;
  } else {
    w = Rcpp::as<arma::mat>(W["r"]);
    if(trans){
      tt = false;
    } else {
      tt = true;
    }
  }
  m = w.n_cols;
  n = u.n_cols;
  for(int i = 0; i < n; i++){
    U = u.col(i);
    U.reshape(m, m);
    U.diag() = 0.5 * U.diag();
    for(int k = 0; k < m; k++){
      for(int r = 0; r < k; r++){
	U(r, k) = 0.0;
      }
    } 
    if(tt){
      a = U * w;
      ans = w.t() * a + a.t() * w;
    } else {
      a = w * U;
      ans = w * a.t() + a * w.t();
    }
    ans.reshape(m * m, 1);
    u.col(i) = ans;
  }
  return(Rcpp::wrap(u));
}



