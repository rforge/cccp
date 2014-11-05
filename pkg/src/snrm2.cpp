/*
 * Norm of a vector in S.
 * snrm2_nls is used for nonlinear, linear and second-order cone constraints.
 * snrm2_p is used for positive semi-definite constraints.
*/

#include "Rcpp.h"
#include "cccp.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(".snrm2_nls")]]
double snrm2_nls(SEXP us) {
  double ans = sqrt(sdot_nls(us, us));
  return ans;
}

// [[Rcpp::export(".snrm2_p")]]
double snrm2_p(SEXP us, SEXP ms) {
  double ans = sqrt(sdot_p(us, us, ms));
  return ans;
}
