#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include "cccp.h"

// [[Rcpp::export(".sslb_nl")]]
SEXP sslb_nl(SEXP us, SEXP lambdas, SEXP invs) {

  Rcpp::NumericMatrix u(us);
  Rcpp::NumericMatrix l(lambdas);
  bool inv = Rcpp::as<bool>(invs);
  int m = u.nrow();

  if(inv){
    for(int i = 0; i < m; i++){
      u(i, 0) *= l(i, 0);
    }
  } else {
    for(int i = 0; i < m; i++){
      u(i, 0) /= l(i, 0);
    }
  }

  return(Rcpp::wrap(u));
}

// [[Rcpp::export(".sslb_s")]]
SEXP sslb_s(SEXP us, SEXP lambdas, SEXP invs) {

  arma::mat u = Rcpp::as<arma::mat>(us);
  arma::mat l = Rcpp::as<arma::mat>(lambdas);
  bool inv = Rcpp::as<bool>(invs);
  double a, cc, lx, u0;
  int m = u.n_rows;

  a = jnrm2(lambdas);
  if(inv){
    lx = dot(l, u) / a;
  } else {
    lx = jdot(lambdas, us) / a;
  }
  u0 = u(0, 0);
  u(0, 0) = lx;
  cc = (lx + u0) / (l(0, 0) / a + 1) / a;
  if(inv == false){
    cc = -cc;
    a = 1 / a;
  }
  for(int i = 1; i < m; i++){
    u(i, 0) = cc * l(i, 0) + u(i, 0);
  }
  u = a * u; 

  return(Rcpp::wrap(u));
}

// [[Rcpp::export(".sslb_p")]]
SEXP sslb_p(SEXP us, SEXP lambdas, SEXP invs, SEXP ms) {

  arma::mat u = Rcpp::as<arma::mat>(us);
  arma::mat l = Rcpp::as<arma::mat>(lambdas);
  arma::vec ld, ls;
  bool inv = Rcpp::as<bool>(invs);
  int m = Rcpp::as<int>(ms);

  u.reshape(m, m);
  l.reshape(m, m);
  ld = l.diag();
  for(int i = 0; i < m; i++){
    ls = sqrt(ld(i)) * sqrt(ld);
    if(inv){
      u.col(i) = u.col(i) % ls;
    } else {
      u.col(i) = u.col(i) / ls;
    }
  }
  u.reshape(m * m, 1);

  return(Rcpp::wrap(u));
}
