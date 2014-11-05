/*
 * Returns x' * J * y, whereby J = [1, 0; 0, -I]
*/

#include "Rcpp.h"

// [[Rcpp::export(".jdot")]]
double jdot(SEXP us, SEXP vs) {

  Rcpp::NumericMatrix u(us);
  Rcpp::NumericMatrix v(vs);
  double ans = 0.0;
  int n = u.nrow();

  ans = u(0, 0) * v(0, 0); 
  for(int i = 1; i < n; i++){
    ans -= u(i, 0) * v(i, 0);
  }

  return(ans);
}
