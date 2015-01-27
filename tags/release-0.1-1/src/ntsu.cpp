/*
 * Update of Nesterov-Todd Scalings
 */

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "cccp.h"

// [[Rcpp::export(".ntsu_n")]]
SEXP ntsu_n(SEXP Ws, SEXP ss, SEXP zs) {

  Rcpp::NumericMatrix s(Rcpp::clone(ss));
  Rcpp::NumericMatrix z(Rcpp::clone(zs));
  int n = s.nrow();
  Rcpp::List W(Rcpp::clone(Ws));
  Rcpp::NumericMatrix dnl = Rcpp::as<Rcpp::NumericMatrix>(W["dnl"]);
  Rcpp::NumericMatrix dnli = Rcpp::as<Rcpp::NumericMatrix>(W["dnli"]);
  Rcpp::S4 lambda = Rcpp::as<Rcpp::S4>(W["lambda"]);
  Rcpp::NumericMatrix l = Rcpp::as<Rcpp::NumericMatrix>(lambda.slot("u"));
  Rcpp::Language NLFS_create("new", "NLFS");
  Rcpp::S4 ans(NLFS_create.eval());
  double ssqrt, zsqrt;

  for(int i = 0; i < n; i++){
    ssqrt = sqrt(s(i, 0));
    zsqrt = sqrt(z(i, 0));

    dnl(i, 0) = dnl(i, 0) * ssqrt / zsqrt;
    dnli(i, 0) = 1.0 / dnl(i, 0);
    l(i, 0) = ssqrt * zsqrt;
  }

  lambda.slot("u") = Rcpp::wrap(l);
  Rcpp::List Wl = Rcpp::List::create(Rcpp::Named("dnl") = dnl,
				     Rcpp::Named("dnli") = dnli,
				     Rcpp::Named("lambda") = lambda);
  ans.slot("W") = Wl;
 
  return(Rcpp::wrap(ans));
}

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

// [[Rcpp::export(".ntsu_s")]]
SEXP ntsu_s(SEXP Ws, SEXP ss, SEXP zs) {

  arma::mat sa = Rcpp::as<arma::mat>(ss);
  arma::mat za = Rcpp::as<arma::mat>(zs);
  int n = za.n_rows;
  Rcpp::List W(Ws);
  arma::mat va = Rcpp::as<arma::mat>(W["v"]);
  double beta = Rcpp::as<double>(W["beta"]);
  Rcpp::S4 lambda = Rcpp::as<Rcpp::S4>(W["lambda"]);
  arma::mat la = Rcpp::as<arma::mat>(lambda.slot("u"));
  Rcpp::Language SOCS_create("new", "SOCS");
  Rcpp::S4 ans(SOCS_create.eval());
  double aa, bb, cc, dd, vs, vz, vq, vu, wk0;

  aa = jnrm2(ss);
  bb = jnrm2(zs);
  sa /=  aa;
  za /=  bb;
  cc = sqrt((1 + arma::dot(sa, za)) / 2.0);
  vs = arma::dot(va, sa);
  vz = va(0, 0) * za(0, 0);
  for(int i = 1; i < n; i++){
    vz -= va(i, 0) * za(i, 0);
  }
  vq = (vs + vz) / 2.0 / cc;
  vu = vs - vz;
  la(0, 0) = cc;
  wk0 = 2 * va(0, 0) * vq - (sa(0, 0) + za(0, 0)) / 2.0 / cc;
  dd = (va(0, 0) * vu - sa(0, 0) / 2.0 + za(0, 0) / 2.0) / (wk0 + 1.0);
  for(int i = 1; i < n; i++){
    la(i, 0) = va(i, 0);
    la(i, 0) *= 2.0 * (-dd * vq + 0.5 * vu); 
    la(i, 0) += 0.5 * (1.0 - dd / cc) * sa(i, 0); 
    la(i, 0) += 0.5 * (1.0 + dd / cc) * za(i, 0); 
  }
  la *= sqrt(aa * bb);
  va *= 2.0 * vq;
  va(0, 0) -= sa(0, 0) / 2.0 / cc;
  for(int i = 1; i < n; i++){
    va(i, 0) += 0.5 / cc * sa(i, 0);
  }
  va += -0.5 / cc * za;
  va(0, 0) += 1.0;
  va *= 1.0 / sqrt(2.0 * va(0, 0));
  beta *= sqrt(aa / bb);

  lambda.slot("u") = Rcpp::wrap(la);
  Rcpp::List Wl = Rcpp::List::create(Rcpp::Named("v") = Rcpp::wrap(va),
				     Rcpp::Named("beta") = beta,
				     Rcpp::Named("lambda") = lambda);
  ans.slot("W") = Wl;
 
  return(Rcpp::wrap(ans));
}

// [[Rcpp::export(".ntsu_p")]]
SEXP ntsu_p(SEXP Ws, SEXP ss, SEXP zs, SEXP ms) {

  arma::mat sa = Rcpp::as<arma::mat>(ss);
  arma::mat za = Rcpp::as<arma::mat>(zs);
  int m = Rcpp::as<int>(ms);
  Rcpp::List W(Ws);
  Rcpp::S4 lambda = Rcpp::as<Rcpp::S4>(W["lambda"]);
  arma::mat ra = Rcpp::as<arma::mat>(W["r"]);
  arma::mat rtia = Rcpp::as<arma::mat>(W["rti"]);
  Rcpp::Language PSDS_create("new", "PSDS");
  Rcpp::S4 ans(PSDS_create.eval());
  arma::mat zts, U, V, DiagL, lu;
  arma::vec l;

  sa.reshape(m, m);
  za.reshape(m, m);
  zts = za.t() * sa;
  svd(U, l, V, zts);
  DiagL = diagmat(1.0 / sqrt(l));
  ra = ra * sa * V * DiagL;
  rtia = rtia * za * U * DiagL;
  lu = diagmat(l);
  lu.reshape(m * m, 1);
  lambda.slot("u") = Rcpp::wrap(lu); 
  Rcpp::List Wl = Rcpp::List::create(Rcpp::Named("r") = Rcpp::wrap(ra),
				     Rcpp::Named("rti") = Rcpp::wrap(rtia),
				     Rcpp::Named("lambda") = lambda);
  ans.slot("W") = Wl; 

  return(Rcpp::wrap(ans));
}

