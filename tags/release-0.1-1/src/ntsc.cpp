/*
 * Computes Nesterov-Todd Scalings
 */

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "cccp.h"

// [[Rcpp::export(".ntsc_n")]]
SEXP ntsc_n(SEXP ss, SEXP zs) {

  Rcpp::NumericMatrix s(ss);
  Rcpp::NumericMatrix z(zs);
  int m = z.nrow();
  Rcpp::Language NLFS_create("new", "NLFS");
  Rcpp::S4 ans(NLFS_create.eval());
  Rcpp::Language NLFV_create("new", "NLFV");
  Rcpp::S4 lambda(NLFV_create.eval());
  Rcpp::NumericMatrix dnl(m, 1), dnli(m, 1), lu(m, 1);

  for(int i = 0; i < m; i++){
    dnl(i, 0) = sqrt(s(i, 0) / z(i, 0));
    dnli(i, 0) = sqrt(z(i, 0) / s(i, 0));
    lu(i, 0) = sqrt(s(i, 0) * z(i, 0));
  }
  lambda.slot("u") = lu; 
  lambda.slot("dims") = m;
  Rcpp::List Wl = Rcpp::List::create(Rcpp::Named("dnl") = dnl,
				     Rcpp::Named("dnli") = dnli,
				     Rcpp::Named("lambda") = lambda);
  ans.slot("W") = Wl; 

  return(Rcpp::wrap(ans));
}

// [[Rcpp::export(".ntsc_l")]]
SEXP ntsc_l(SEXP ss, SEXP zs) {

  Rcpp::NumericMatrix s(ss);
  Rcpp::NumericMatrix z(zs);
  int m = z.nrow();
  Rcpp::Language NNOS_create("new", "NNOS");
  Rcpp::S4 ans(NNOS_create.eval());
  Rcpp::Language NNOV_create("new", "NNOV");
  Rcpp::S4 lambda(NNOV_create.eval());
  Rcpp::NumericMatrix d(m, 1), di(m, 1), lu(m, 1);

  for(int i = 0; i < m; i++){
    d(i, 0) = sqrt(s(i, 0) / z(i, 0));
    di(i, 0) = sqrt(z(i, 0) / s(i, 0));
    lu(i, 0) = sqrt(s(i, 0) * z(i, 0));
  }
  lambda.slot("u") = lu; 
  lambda.slot("dims") = m;
  Rcpp::List Wl = Rcpp::List::create(Rcpp::Named("d") = d,
				     Rcpp::Named("di") = di,
				     Rcpp::Named("lambda") = lambda);
  ans.slot("W") = Wl; 

  return(Rcpp::wrap(ans));
}

// [[Rcpp::export(".ntsc_s")]]
SEXP ntsc_s(SEXP ss, SEXP zs) {

  Rcpp::NumericMatrix s(ss);
  int m = s.nrow();
  Rcpp::NumericMatrix z(zs);
  Rcpp::Language SOCS_create("new", "SOCS");
  Rcpp::S4 ans(SOCS_create.eval());
  Rcpp::Language SOCV_create("new", "SOCV");
  Rcpp::S4 lambda(SOCV_create.eval());
  arma::mat v(m, 1), lu(m, 1);
  double aa, bb, cc, dd, beta, szdot;

  aa = jnrm2(ss);
  bb = jnrm2(zs);
  beta = sqrt(aa / bb);
  szdot = sdot_nls(ss, zs);
  cc = sqrt((szdot / aa / bb + 1.0) / 2.0);
  v = -z / bb;
  v(0, 0) = -v(0, 0);
  for(int i = 0; i < m; i++){
    v(i, 0) += s(i, 0) / aa;
    v(i, 0) *= (1.0 / 2.0 / cc);
  }
  v(0, 0) = v(0, 0) + 1.0;
  v *= 1.0 / sqrt(2.0 * v(0, 0));
  lu(0, 0) = cc;
  dd = 2 * cc + s(0, 0) / aa + z(0, 0) / bb;
  for(int i = 1; i < m; i++){
    lu(i, 0) = s(i, 0);
    lu(i, 0) = (cc + z(0, 0) / bb) / dd / aa * lu(i, 0);
    lu(i, 0) = (cc + s(0, 0) / aa) / dd / bb * z(i, 0) + lu(i, 0); 
  }
  lambda.slot("u") = Rcpp::wrap(sqrt(aa * bb) * lu); 
  lambda.slot("dims") = m;
  Rcpp::List Wl = Rcpp::List::create(Rcpp::Named("v") = Rcpp::wrap(v),
				     Rcpp::Named("beta") = Rcpp::wrap(beta),
				     Rcpp::Named("lambda") = lambda);
  ans.slot("W") = Wl; 

  return(Rcpp::wrap(ans));
}

// [[Rcpp::export(".ntsc_p")]]
SEXP ntsc_p(SEXP ss, SEXP zs, SEXP ms) {

  arma::mat s = Rcpp::as<arma::mat>(ss);
  arma::mat z = Rcpp::as<arma::mat>(zs);
  int m = Rcpp::as<int>(ms);
  Rcpp::Language PSDS_create("new", "PSDS");
  Rcpp::S4 ans(PSDS_create.eval());
  Rcpp::Language PSDV_create("new", "PSDV");
  Rcpp::S4 lambda(PSDV_create.eval());
  arma::mat sc, zc, szc, U, V;
  arma::vec l;
  arma::mat r, rti, lu;

  s.reshape(m, m);
  z.reshape(m, m);
  arma::chol(sc, s);
  arma::chol(zc, z);
  szc = zc * sc.t();
  svd(U, l, V, szc);
  r = zc.i() * U * diagmat(sqrt(l));
  rti = zc.t() * U * diagmat(1.0 / sqrt(l));
  lu = diagmat(l);
  lu.reshape(m * m, 1);

  lambda.slot("u") = Rcpp::wrap(lu); 
  lambda.slot("dims") = m;
  Rcpp::List Wl = Rcpp::List::create(Rcpp::Named("r") = Rcpp::wrap(r),
				     Rcpp::Named("rti") = Rcpp::wrap(rti),
				     Rcpp::Named("lambda") = lambda);
  ans.slot("W") = Wl; 

  return(Rcpp::wrap(ans));
}

