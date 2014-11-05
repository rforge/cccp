// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// jdot
double jdot(SEXP us, SEXP vs);
RcppExport SEXP cccp_jdot(SEXP usSEXP, SEXP vsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type vs(vsSEXP );
        double __result = jdot(us, vs);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// jnrm2
double jnrm2(SEXP us);
RcppExport SEXP cccp_jnrm2(SEXP usSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        double __result = jnrm2(us);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// sdot_nls
double sdot_nls(SEXP us, SEXP vs);
RcppExport SEXP cccp_sdot_nls(SEXP usSEXP, SEXP vsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type vs(vsSEXP );
        double __result = sdot_nls(us, vs);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// sdot_p
double sdot_p(SEXP us, SEXP vs, SEXP ms);
RcppExport SEXP cccp_sdot_p(SEXP usSEXP, SEXP vsSEXP, SEXP msSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type vs(vsSEXP );
        Rcpp::traits::input_parameter< SEXP >::type ms(msSEXP );
        double __result = sdot_p(us, vs, ms);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// sinv_nl
SEXP sinv_nl(SEXP us, SEXP vs);
RcppExport SEXP cccp_sinv_nl(SEXP usSEXP, SEXP vsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type vs(vsSEXP );
        SEXP __result = sinv_nl(us, vs);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// sinv_s
SEXP sinv_s(SEXP us, SEXP vs);
RcppExport SEXP cccp_sinv_s(SEXP usSEXP, SEXP vsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type vs(vsSEXP );
        SEXP __result = sinv_s(us, vs);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// sinv_p
SEXP sinv_p(SEXP us, SEXP vs, SEXP ms);
RcppExport SEXP cccp_sinv_p(SEXP usSEXP, SEXP vsSEXP, SEXP msSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type vs(vsSEXP );
        Rcpp::traits::input_parameter< SEXP >::type ms(msSEXP );
        SEXP __result = sinv_p(us, vs, ms);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// smsa_nl
SEXP smsa_nl(SEXP us, SEXP alphas, SEXP inits);
RcppExport SEXP cccp_smsa_nl(SEXP usSEXP, SEXP alphasSEXP, SEXP initsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type alphas(alphasSEXP );
        Rcpp::traits::input_parameter< SEXP >::type inits(initsSEXP );
        SEXP __result = smsa_nl(us, alphas, inits);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// smsa_s
SEXP smsa_s(SEXP us, SEXP alphas, SEXP inits);
RcppExport SEXP cccp_smsa_s(SEXP usSEXP, SEXP alphasSEXP, SEXP initsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type alphas(alphasSEXP );
        Rcpp::traits::input_parameter< SEXP >::type inits(initsSEXP );
        SEXP __result = smsa_s(us, alphas, inits);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// smsa1_p
SEXP smsa1_p(SEXP us, SEXP alphas, SEXP ms);
RcppExport SEXP cccp_smsa1_p(SEXP usSEXP, SEXP alphasSEXP, SEXP msSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type alphas(alphasSEXP );
        Rcpp::traits::input_parameter< SEXP >::type ms(msSEXP );
        SEXP __result = smsa1_p(us, alphas, ms);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// smsa2_p
SEXP smsa2_p(SEXP us, SEXP alphas, SEXP sigmas, SEXP lambdas, SEXP ms);
RcppExport SEXP cccp_smsa2_p(SEXP usSEXP, SEXP alphasSEXP, SEXP sigmasSEXP, SEXP lambdasSEXP, SEXP msSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type alphas(alphasSEXP );
        Rcpp::traits::input_parameter< SEXP >::type sigmas(sigmasSEXP );
        Rcpp::traits::input_parameter< SEXP >::type lambdas(lambdasSEXP );
        Rcpp::traits::input_parameter< SEXP >::type ms(msSEXP );
        SEXP __result = smsa2_p(us, alphas, sigmas, lambdas, ms);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// smss_nl
SEXP smss_nl(SEXP us);
RcppExport SEXP cccp_smss_nl(SEXP usSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        SEXP __result = smss_nl(us);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// smss_s
SEXP smss_s(SEXP us);
RcppExport SEXP cccp_smss_s(SEXP usSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        SEXP __result = smss_s(us);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// smss_p
SEXP smss_p(SEXP us, SEXP ms);
RcppExport SEXP cccp_smss_p(SEXP usSEXP, SEXP msSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type ms(msSEXP );
        SEXP __result = smss_p(us, ms);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// snrm2_nls
double snrm2_nls(SEXP us);
RcppExport SEXP cccp_snrm2_nls(SEXP usSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        double __result = snrm2_nls(us);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// snrm2_p
double snrm2_p(SEXP us, SEXP ms);
RcppExport SEXP cccp_snrm2_p(SEXP usSEXP, SEXP msSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type ms(msSEXP );
        double __result = snrm2_p(us, ms);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// sone_nl
SEXP sone_nl(SEXP us);
RcppExport SEXP cccp_sone_nl(SEXP usSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        SEXP __result = sone_nl(us);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// sone_s
SEXP sone_s(SEXP us);
RcppExport SEXP cccp_sone_s(SEXP usSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        SEXP __result = sone_s(us);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// sone_p
SEXP sone_p(SEXP ms);
RcppExport SEXP cccp_sone_p(SEXP msSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type ms(msSEXP );
        SEXP __result = sone_p(ms);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// sprd_nl
SEXP sprd_nl(SEXP us, SEXP vs);
RcppExport SEXP cccp_sprd_nl(SEXP usSEXP, SEXP vsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type vs(vsSEXP );
        SEXP __result = sprd_nl(us, vs);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// sprd_s
SEXP sprd_s(SEXP us, SEXP vs);
RcppExport SEXP cccp_sprd_s(SEXP usSEXP, SEXP vsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type vs(vsSEXP );
        SEXP __result = sprd_s(us, vs);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// sprd_p
SEXP sprd_p(SEXP us, SEXP vs, SEXP ms);
RcppExport SEXP cccp_sprd_p(SEXP usSEXP, SEXP vsSEXP, SEXP msSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type vs(vsSEXP );
        Rcpp::traits::input_parameter< SEXP >::type ms(msSEXP );
        SEXP __result = sprd_p(us, vs, ms);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// ssnt_l
SEXP ssnt_l(SEXP us, SEXP Ws, SEXP invs);
RcppExport SEXP cccp_ssnt_l(SEXP usSEXP, SEXP WsSEXP, SEXP invsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type Ws(WsSEXP );
        Rcpp::traits::input_parameter< SEXP >::type invs(invsSEXP );
        SEXP __result = ssnt_l(us, Ws, invs);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// ssnt_s
SEXP ssnt_s(SEXP us, SEXP Ws, SEXP invs);
RcppExport SEXP cccp_ssnt_s(SEXP usSEXP, SEXP WsSEXP, SEXP invsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type us(usSEXP );
        Rcpp::traits::input_parameter< SEXP >::type Ws(WsSEXP );
        Rcpp::traits::input_parameter< SEXP >::type invs(invsSEXP );
        SEXP __result = ssnt_s(us, Ws, invs);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
