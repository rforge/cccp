/*
 * Returns sqrt(x' * J * y), whereby J = [1, 0; 0, -I]
*/

#include "Rcpp.h"
#include "cccp.h"

// [[Rcpp::export(".jnrm2")]]
double jnrm2(SEXP us) {
  double ans = sqrt(jdot(us, us));
  return ans;
}
