/* 
 *
 * Header file for package cccp
 *
*/

#ifndef CCCP_H
#define CCCP_H

/*
 * Functions for manipulating vectors in S
*/

double sdot_nls(SEXP us, SEXP vs);
double sdot_p(SEXP us, SEXP vs, SEXP ms);
double snrm2_nls(SEXP us);
double snrm2_p(SEXP us, SEXP ms);
double jdot(SEXP us, SEXP vs);
double jnrm2(SEXP us);
SEXP sprd_nl(SEXP us, SEXP vs);
SEXP sprd_s(SEXP us, SEXP vs);
SEXP sprd_p(SEXP us, SEXP vs, SEXP ms);
SEXP sinv_nl(SEXP us, SEXP vs);
SEXP sinv_s(SEXP us, SEXP vs);
SEXP sinv_p(SEXP us, SEXP vs, SEXP ms);
SEXP smss_nl(SEXP us);
SEXP smss_s(SEXP us);
SEXP smss_p(SEXP us, SEXP ms);
SEXP sone_nl(SEXP us);
SEXP sone_s(SEXP us);
SEXP sone_p(SEXP ms);
SEXP smsa_nl(SEXP us, SEXP alphas, SEXP inits);
SEXP smsa_s(SEXP us, SEXP alphas, SEXP inits);
SEXP smsa1_p(SEXP us, SEXP alphas, SEXP ms);
SEXP smsa2_p(SEXP us, SEXP alphas, SEXP sigmas, SEXP lambdas, SEXP ms);
SEXP ssnt_n(SEXP us, SEXP Ws, SEXP invs);
SEXP ssnt_l(SEXP us, SEXP Ws, SEXP invs);
SEXP ssnt_s(SEXP us, SEXP Ws, SEXP invs);
SEXP ssnt_p(SEXP us, SEXP Ws, SEXP transs, SEXP invs); 
SEXP sslb_nl(SEXP us, SEXP lambdas, SEXP invs);
SEXP sslb_s(SEXP us, SEXP lambdas, SEXP invs);
SEXP sslb_p(SEXP us, SEXP lambdas, SEXP invs, SEXP ms);
SEXP ntsc_n(SEXP ss, SEXP zs);
SEXP ntsc_l(SEXP ss, SEXP zs);
SEXP ntsc_s(SEXP ss, SEXP zs);
SEXP ntsc_p(SEXP ss, SEXP zs, SEXP ms);

#endif
