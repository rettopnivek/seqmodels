// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// pdiff
Rcpp::NumericVector pdiff(Rcpp::NumericVector rt, Rcpp::NumericVector ch, Rcpp::NumericVector alpha, Rcpp::NumericVector theta, Rcpp::NumericVector xi, Rcpp::NumericVector tau, Rcpp::NumericVector sigma, double eps, int parYes);
RcppExport SEXP seqmodels_pdiff(SEXP rtSEXP, SEXP chSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP xiSEXP, SEXP tauSEXP, SEXP sigmaSEXP, SEXP epsSEXP, SEXP parYesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ch(chSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type parYes(parYesSEXP);
    __result = Rcpp::wrap(pdiff(rt, ch, alpha, theta, xi, tau, sigma, eps, parYes));
    return __result;
END_RCPP
}
// ddiff
Rcpp::NumericVector ddiff(Rcpp::NumericVector rt, Rcpp::NumericVector ch, Rcpp::NumericVector alpha, Rcpp::NumericVector theta, Rcpp::NumericVector xi, Rcpp::NumericVector tau, Rcpp::NumericVector eta, Rcpp::NumericVector stheta, Rcpp::NumericVector stau, Rcpp::NumericVector sigma, double eps, int ln, int parYes);
RcppExport SEXP seqmodels_ddiff(SEXP rtSEXP, SEXP chSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP xiSEXP, SEXP tauSEXP, SEXP etaSEXP, SEXP sthetaSEXP, SEXP stauSEXP, SEXP sigmaSEXP, SEXP epsSEXP, SEXP lnSEXP, SEXP parYesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ch(chSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type stheta(sthetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type stau(stauSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type ln(lnSEXP);
    Rcpp::traits::input_parameter< int >::type parYes(parYesSEXP);
    __result = Rcpp::wrap(ddiff(rt, ch, alpha, theta, xi, tau, eta, stheta, stau, sigma, eps, ln, parYes));
    return __result;
END_RCPP
}
// rinvgauss
Rcpp::NumericVector rinvgauss(int N, Rcpp::NumericVector kappa, Rcpp::NumericVector xi, Rcpp::NumericVector sigma);
RcppExport SEXP seqmodels_rinvgauss(SEXP NSEXP, SEXP kappaSEXP, SEXP xiSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(rinvgauss(N, kappa, xi, sigma));
    return __result;
END_RCPP
}
// dinvgauss
Rcpp::NumericVector dinvgauss(Rcpp::NumericVector t, Rcpp::NumericVector kappa, Rcpp::NumericVector xi, Rcpp::NumericVector sigma, int ln);
RcppExport SEXP seqmodels_dinvgauss(SEXP tSEXP, SEXP kappaSEXP, SEXP xiSEXP, SEXP sigmaSEXP, SEXP lnSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type ln(lnSEXP);
    __result = Rcpp::wrap(dinvgauss(t, kappa, xi, sigma, ln));
    return __result;
END_RCPP
}
// pinvgauss
Rcpp::NumericVector pinvgauss(Rcpp::NumericVector t, Rcpp::NumericVector kappa, Rcpp::NumericVector xi, Rcpp::NumericVector sigma);
RcppExport SEXP seqmodels_pinvgauss(SEXP tSEXP, SEXP kappaSEXP, SEXP xiSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(pinvgauss(t, kappa, xi, sigma));
    return __result;
END_RCPP
}
// rwaldrace
Rcpp::NumericMatrix rwaldrace(int N, Rcpp::NumericVector k1, Rcpp::NumericVector xi1, Rcpp::NumericVector tau1, Rcpp::NumericVector k0, Rcpp::NumericVector xi0, Rcpp::NumericVector tau0, Rcpp::NumericVector s1, Rcpp::NumericVector s0);
RcppExport SEXP seqmodels_rwaldrace(SEXP NSEXP, SEXP k1SEXP, SEXP xi1SEXP, SEXP tau1SEXP, SEXP k0SEXP, SEXP xi0SEXP, SEXP tau0SEXP, SEXP s1SEXP, SEXP s0SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi1(xi1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau1(tau1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi0(xi0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau0(tau0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type s0(s0SEXP);
    __result = Rcpp::wrap(rwaldrace(N, k1, xi1, tau1, k0, xi0, tau0, s1, s0));
    return __result;
END_RCPP
}
// dwaldrace
Rcpp::NumericVector dwaldrace(Rcpp::NumericVector rt, Rcpp::NumericVector ch, Rcpp::NumericVector k1, Rcpp::NumericVector xi1, Rcpp::NumericVector tau1, Rcpp::NumericVector k0, Rcpp::NumericVector xi0, Rcpp::NumericVector tau0, Rcpp::NumericVector s1, Rcpp::NumericVector s0, int ln);
RcppExport SEXP seqmodels_dwaldrace(SEXP rtSEXP, SEXP chSEXP, SEXP k1SEXP, SEXP xi1SEXP, SEXP tau1SEXP, SEXP k0SEXP, SEXP xi0SEXP, SEXP tau0SEXP, SEXP s1SEXP, SEXP s0SEXP, SEXP lnSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ch(chSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi1(xi1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau1(tau1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi0(xi0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau0(tau0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type s0(s0SEXP);
    Rcpp::traits::input_parameter< int >::type ln(lnSEXP);
    __result = Rcpp::wrap(dwaldrace(rt, ch, k1, xi1, tau1, k0, xi0, tau0, s1, s0, ln));
    return __result;
END_RCPP
}
// pwaldrace
Rcpp::NumericVector pwaldrace(Rcpp::NumericVector rt, Rcpp::NumericVector ch, Rcpp::NumericVector k1, Rcpp::NumericVector xi1, Rcpp::NumericVector tau1, Rcpp::NumericVector k0, Rcpp::NumericVector xi0, Rcpp::NumericVector tau0, Rcpp::NumericVector s1, Rcpp::NumericVector s0);
RcppExport SEXP seqmodels_pwaldrace(SEXP rtSEXP, SEXP chSEXP, SEXP k1SEXP, SEXP xi1SEXP, SEXP tau1SEXP, SEXP k0SEXP, SEXP xi0SEXP, SEXP tau0SEXP, SEXP s1SEXP, SEXP s0SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ch(chSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi1(xi1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau1(tau1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi0(xi0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau0(tau0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type s0(s0SEXP);
    __result = Rcpp::wrap(pwaldrace(rt, ch, k1, xi1, tau1, k0, xi0, tau0, s1, s0));
    return __result;
END_RCPP
}
// pwaldrace2
Rcpp::NumericVector pwaldrace2(Rcpp::NumericVector rt, Rcpp::NumericVector ch, Rcpp::NumericVector k1, Rcpp::NumericVector xi1, Rcpp::NumericVector tau1, Rcpp::NumericVector k0, Rcpp::NumericVector xi0, Rcpp::NumericVector tau0, Rcpp::NumericVector s1, Rcpp::NumericVector s0, int parYes);
RcppExport SEXP seqmodels_pwaldrace2(SEXP rtSEXP, SEXP chSEXP, SEXP k1SEXP, SEXP xi1SEXP, SEXP tau1SEXP, SEXP k0SEXP, SEXP xi0SEXP, SEXP tau0SEXP, SEXP s1SEXP, SEXP s0SEXP, SEXP parYesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ch(chSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi1(xi1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau1(tau1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi0(xi0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau0(tau0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type s0(s0SEXP);
    Rcpp::traits::input_parameter< int >::type parYes(parYesSEXP);
    __result = Rcpp::wrap(pwaldrace2(rt, ch, k1, xi1, tau1, k0, xi0, tau0, s1, s0, parYes));
    return __result;
END_RCPP
}
