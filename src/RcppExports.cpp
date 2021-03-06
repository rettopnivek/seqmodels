// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// remg
Rcpp::NumericVector remg(int n, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector lambda);
RcppExport SEXP _seqmodels_remg(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(remg(n, mu, sigma, lambda));
    return rcpp_result_gen;
END_RCPP
}
// demg
Rcpp::NumericVector demg(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector lambda, bool ln);
RcppExport SEXP _seqmodels_demg(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP lambdaSEXP, SEXP lnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< bool >::type ln(lnSEXP);
    rcpp_result_gen = Rcpp::wrap(demg(x, mu, sigma, lambda, ln));
    return rcpp_result_gen;
END_RCPP
}
// pemg
Rcpp::NumericVector pemg(Rcpp::NumericVector q, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector lambda, bool ln, bool lower_tail);
RcppExport SEXP _seqmodels_pemg(SEXP qSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP lambdaSEXP, SEXP lnSEXP, SEXP lower_tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< bool >::type ln(lnSEXP);
    Rcpp::traits::input_parameter< bool >::type lower_tail(lower_tailSEXP);
    rcpp_result_gen = Rcpp::wrap(pemg(q, mu, sigma, lambda, ln, lower_tail));
    return rcpp_result_gen;
END_RCPP
}
// qemg
Rcpp::NumericVector qemg(Rcpp::NumericVector p, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector lambda, Rcpp::NumericVector bounds, double em_stop, double err);
RcppExport SEXP _seqmodels_qemg(SEXP pSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP lambdaSEXP, SEXP boundsSEXP, SEXP em_stopSEXP, SEXP errSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< double >::type em_stop(em_stopSEXP);
    Rcpp::traits::input_parameter< double >::type err(errSEXP);
    rcpp_result_gen = Rcpp::wrap(qemg(p, mu, sigma, lambda, bounds, em_stop, err));
    return rcpp_result_gen;
END_RCPP
}
// memg
Rcpp::DataFrame memg(Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector lambda);
RcppExport SEXP _seqmodels_memg(SEXP muSEXP, SEXP sigmaSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(memg(mu, sigma, lambda));
    return rcpp_result_gen;
END_RCPP
}
// dlevy
Rcpp::NumericVector dlevy(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, bool ln);
RcppExport SEXP _seqmodels_dlevy(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP lnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type ln(lnSEXP);
    rcpp_result_gen = Rcpp::wrap(dlevy(x, mu, sigma, ln));
    return rcpp_result_gen;
END_RCPP
}
// plevy
Rcpp::NumericVector plevy(Rcpp::NumericVector q, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, bool lower_tail, bool ln);
RcppExport SEXP _seqmodels_plevy(SEXP qSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP lower_tailSEXP, SEXP lnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type lower_tail(lower_tailSEXP);
    Rcpp::traits::input_parameter< bool >::type ln(lnSEXP);
    rcpp_result_gen = Rcpp::wrap(plevy(q, mu, sigma, lower_tail, ln));
    return rcpp_result_gen;
END_RCPP
}
// qlevy
Rcpp::NumericVector qlevy(Rcpp::NumericVector p, Rcpp::NumericVector mu, Rcpp::NumericVector sigma);
RcppExport SEXP _seqmodels_qlevy(SEXP pSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(qlevy(p, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// rlevy
Rcpp::NumericVector rlevy(int n, Rcpp::NumericVector mu, Rcpp::NumericVector sigma);
RcppExport SEXP _seqmodels_rlevy(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(rlevy(n, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// rinvgauss
Rcpp::NumericVector rinvgauss(int n, Rcpp::NumericVector kappa, Rcpp::NumericVector xi, Rcpp::NumericVector tau, Rcpp::NumericVector sigma);
RcppExport SEXP _seqmodels_rinvgauss(SEXP nSEXP, SEXP kappaSEXP, SEXP xiSEXP, SEXP tauSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(rinvgauss(n, kappa, xi, tau, sigma));
    return rcpp_result_gen;
END_RCPP
}
// dinvgauss
Rcpp::NumericVector dinvgauss(Rcpp::NumericVector t, Rcpp::NumericVector kappa, Rcpp::NumericVector xi, Rcpp::NumericVector tau, Rcpp::NumericVector sigma, bool ln);
RcppExport SEXP _seqmodels_dinvgauss(SEXP tSEXP, SEXP kappaSEXP, SEXP xiSEXP, SEXP tauSEXP, SEXP sigmaSEXP, SEXP lnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type ln(lnSEXP);
    rcpp_result_gen = Rcpp::wrap(dinvgauss(t, kappa, xi, tau, sigma, ln));
    return rcpp_result_gen;
END_RCPP
}
// pinvgauss
Rcpp::NumericVector pinvgauss(Rcpp::NumericVector t, Rcpp::NumericVector kappa, Rcpp::NumericVector xi, Rcpp::NumericVector tau, Rcpp::NumericVector sigma, bool ln, bool lower_tail);
RcppExport SEXP _seqmodels_pinvgauss(SEXP tSEXP, SEXP kappaSEXP, SEXP xiSEXP, SEXP tauSEXP, SEXP sigmaSEXP, SEXP lnSEXP, SEXP lower_tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type ln(lnSEXP);
    Rcpp::traits::input_parameter< bool >::type lower_tail(lower_tailSEXP);
    rcpp_result_gen = Rcpp::wrap(pinvgauss(t, kappa, xi, tau, sigma, ln, lower_tail));
    return rcpp_result_gen;
END_RCPP
}
// qinvgauss
Rcpp::NumericVector qinvgauss(Rcpp::NumericVector p, Rcpp::NumericVector kappa, Rcpp::NumericVector xi, Rcpp::NumericVector tau, Rcpp::NumericVector sigma, double bounds, double em_stop, double err);
RcppExport SEXP _seqmodels_qinvgauss(SEXP pSEXP, SEXP kappaSEXP, SEXP xiSEXP, SEXP tauSEXP, SEXP sigmaSEXP, SEXP boundsSEXP, SEXP em_stopSEXP, SEXP errSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< double >::type em_stop(em_stopSEXP);
    Rcpp::traits::input_parameter< double >::type err(errSEXP);
    rcpp_result_gen = Rcpp::wrap(qinvgauss(p, kappa, xi, tau, sigma, bounds, em_stop, err));
    return rcpp_result_gen;
END_RCPP
}
// minvgauss
Rcpp::DataFrame minvgauss(Rcpp::NumericVector kappa, Rcpp::NumericVector xi, Rcpp::NumericVector tau, Rcpp::NumericVector sigma);
RcppExport SEXP _seqmodels_minvgauss(SEXP kappaSEXP, SEXP xiSEXP, SEXP tauSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(minvgauss(kappa, xi, tau, sigma));
    return rcpp_result_gen;
END_RCPP
}
// dwiener
Rcpp::NumericVector dwiener(Rcpp::NumericVector rt, Rcpp::NumericVector ch, Rcpp::NumericVector alpha, Rcpp::NumericVector theta, Rcpp::NumericVector xi, Rcpp::NumericVector tau, Rcpp::NumericVector sigma, bool ln, bool joint, double eps, bool parYes);
RcppExport SEXP _seqmodels_dwiener(SEXP rtSEXP, SEXP chSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP xiSEXP, SEXP tauSEXP, SEXP sigmaSEXP, SEXP lnSEXP, SEXP jointSEXP, SEXP epsSEXP, SEXP parYesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ch(chSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type ln(lnSEXP);
    Rcpp::traits::input_parameter< bool >::type joint(jointSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< bool >::type parYes(parYesSEXP);
    rcpp_result_gen = Rcpp::wrap(dwiener(rt, ch, alpha, theta, xi, tau, sigma, ln, joint, eps, parYes));
    return rcpp_result_gen;
END_RCPP
}
// pwiener
Rcpp::NumericVector pwiener(Rcpp::NumericVector rt, Rcpp::NumericVector ch, Rcpp::NumericVector alpha, Rcpp::NumericVector theta, Rcpp::NumericVector xi, Rcpp::NumericVector tau, Rcpp::NumericVector sigma, bool ln, bool joint, bool lower_tail, double eps, bool parYes);
RcppExport SEXP _seqmodels_pwiener(SEXP rtSEXP, SEXP chSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP xiSEXP, SEXP tauSEXP, SEXP sigmaSEXP, SEXP lnSEXP, SEXP jointSEXP, SEXP lower_tailSEXP, SEXP epsSEXP, SEXP parYesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ch(chSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type ln(lnSEXP);
    Rcpp::traits::input_parameter< bool >::type joint(jointSEXP);
    Rcpp::traits::input_parameter< bool >::type lower_tail(lower_tailSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< bool >::type parYes(parYesSEXP);
    rcpp_result_gen = Rcpp::wrap(pwiener(rt, ch, alpha, theta, xi, tau, sigma, ln, joint, lower_tail, eps, parYes));
    return rcpp_result_gen;
END_RCPP
}
// qwiener
Rcpp::NumericVector qwiener(Rcpp::NumericVector p, Rcpp::NumericVector ch, Rcpp::NumericVector alpha, Rcpp::NumericVector theta, Rcpp::NumericVector xi, Rcpp::NumericVector tau, Rcpp::NumericVector sigma, bool joint, double eps, double bounds, double em_stop, double err, bool parYes);
RcppExport SEXP _seqmodels_qwiener(SEXP pSEXP, SEXP chSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP xiSEXP, SEXP tauSEXP, SEXP sigmaSEXP, SEXP jointSEXP, SEXP epsSEXP, SEXP boundsSEXP, SEXP em_stopSEXP, SEXP errSEXP, SEXP parYesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ch(chSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type joint(jointSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< double >::type em_stop(em_stopSEXP);
    Rcpp::traits::input_parameter< double >::type err(errSEXP);
    Rcpp::traits::input_parameter< bool >::type parYes(parYesSEXP);
    rcpp_result_gen = Rcpp::wrap(qwiener(p, ch, alpha, theta, xi, tau, sigma, joint, eps, bounds, em_stop, err, parYes));
    return rcpp_result_gen;
END_RCPP
}
// rwiener
Rcpp::DataFrame rwiener(int n, Rcpp::NumericVector alpha, Rcpp::NumericVector theta, Rcpp::NumericVector xi, Rcpp::NumericVector tau, Rcpp::NumericVector sigma, double eps, double bounds, double em_stop, double err, bool parYes);
RcppExport SEXP _seqmodels_rwiener(SEXP nSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP xiSEXP, SEXP tauSEXP, SEXP sigmaSEXP, SEXP epsSEXP, SEXP boundsSEXP, SEXP em_stopSEXP, SEXP errSEXP, SEXP parYesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< double >::type em_stop(em_stopSEXP);
    Rcpp::traits::input_parameter< double >::type err(errSEXP);
    Rcpp::traits::input_parameter< bool >::type parYes(parYesSEXP);
    rcpp_result_gen = Rcpp::wrap(rwiener(n, alpha, theta, xi, tau, sigma, eps, bounds, em_stop, err, parYes));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_seqmodels_remg", (DL_FUNC) &_seqmodels_remg, 4},
    {"_seqmodels_demg", (DL_FUNC) &_seqmodels_demg, 5},
    {"_seqmodels_pemg", (DL_FUNC) &_seqmodels_pemg, 6},
    {"_seqmodels_qemg", (DL_FUNC) &_seqmodels_qemg, 7},
    {"_seqmodels_memg", (DL_FUNC) &_seqmodels_memg, 3},
    {"_seqmodels_dlevy", (DL_FUNC) &_seqmodels_dlevy, 4},
    {"_seqmodels_plevy", (DL_FUNC) &_seqmodels_plevy, 5},
    {"_seqmodels_qlevy", (DL_FUNC) &_seqmodels_qlevy, 3},
    {"_seqmodels_rlevy", (DL_FUNC) &_seqmodels_rlevy, 3},
    {"_seqmodels_rinvgauss", (DL_FUNC) &_seqmodels_rinvgauss, 5},
    {"_seqmodels_dinvgauss", (DL_FUNC) &_seqmodels_dinvgauss, 6},
    {"_seqmodels_pinvgauss", (DL_FUNC) &_seqmodels_pinvgauss, 7},
    {"_seqmodels_qinvgauss", (DL_FUNC) &_seqmodels_qinvgauss, 8},
    {"_seqmodels_minvgauss", (DL_FUNC) &_seqmodels_minvgauss, 4},
    {"_seqmodels_dwiener", (DL_FUNC) &_seqmodels_dwiener, 11},
    {"_seqmodels_pwiener", (DL_FUNC) &_seqmodels_pwiener, 12},
    {"_seqmodels_qwiener", (DL_FUNC) &_seqmodels_qwiener, 13},
    {"_seqmodels_rwiener", (DL_FUNC) &_seqmodels_rwiener, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_seqmodels(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
