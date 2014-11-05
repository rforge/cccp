#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export(".smss_nl")]]
SEXP smss_nl(SEXP us) {
  Rcpp::NumericVector u(us);
  return(Rcpp::List::create(Rcpp::Named("ms") = -min(u),
			    Rcpp::Named("evd") = R_NilValue));
}

// [[Rcpp::export(".smss_s")]]
SEXP smss_s(SEXP us) {
  Rcpp::NumericMatrix u(us);
  int n = u.nrow();
  double ms = 0.0;

  for(int i = 1; i < n; i++){
    ms += u(i, 0) * u(i, 0);
  }
  ms = sqrt(ms);
  ms = ms - u(0, 0);

  return(Rcpp::List::create(Rcpp::Named("ms") = ms,
			    Rcpp::Named("evd") = R_NilValue));
}

// [[Rcpp::export(".smss_p")]]
SEXP smss_p(SEXP us, SEXP ms) {

  int m = Rcpp::as<int>(ms);
  arma::mat U = Rcpp::as<arma::mat>(us);
  U.reshape(m, m);
  arma::vec eval;
  arma::mat evec;
  arma::eig_sym(eval, evec, U);

  Rcpp::List evd;
  evd = Rcpp::List::create(Rcpp::Named("values") = Rcpp::wrap(eval),
			   Rcpp::Named("vectors") = Rcpp::wrap(evec));


  return(Rcpp::List::create(Rcpp::Named("ms") = Rcpp::wrap(-eval[0]),
			    Rcpp::Named("evd") = evd));
}
