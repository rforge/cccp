/*
 * Update of Nesterov-Todd Scalings
 */

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(".ntsu_l")]]
SEXP ntsu_l(SEXP Ws, SEXP ss, SEXP zs) {

  Rcpp::NumericMatrix s(Rcpp::clone(ss));
  Rcpp::NumericMatrix z(Rcpp::clone(zs));
  int n = s.nrow();
  Rcpp::List W(Rcpp::clone(Ws));
  Rcpp::NumericMatrix d = Rcpp::as<Rcpp::NumericMatrix>(W["d"]);
  Rcpp::NumericMatrix di = Rcpp::as<Rcpp::NumericMatrix>(W["di"]);
  Rcpp::S4 lambda = Rcpp::as<Rcpp::S4>(W["lambda"]);
  Rcpp::NumericMatrix l = Rcpp::as<Rcpp::NumericMatrix>(lambda.slot("u"));
  Rcpp::Language NNOS_create("new", "NNOS");
  Rcpp::S4 ans(NNOS_create.eval());
  double ssqrt, zsqrt;

  for(int i = 0; i < n; i++){
    ssqrt = sqrt(s(i, 0));
    zsqrt = sqrt(z(i, 0));

    d(i, 0) = d(i, 0) * ssqrt / zsqrt;
    di(i, 0) = 1.0 / d(i, 0);
    l(i, 0) = ssqrt * zsqrt;
  }

  lambda.slot("u") = Rcpp::wrap(l);
  Rcpp::List Wl = Rcpp::List::create(Rcpp::Named("d") = d,
				     Rcpp::Named("di") = di,
				     Rcpp::Named("lambda") = lambda);
  ans.slot("W") = Wl;
 
  return(Rcpp::wrap(ans));
}
