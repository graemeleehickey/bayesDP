// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sigma2marginal
SEXP sigma2marginal(int n, const arma::vec& grid, const arma::mat& XtX, const arma::mat& SigmaBetaInv, const arma::mat& Xstar, const arma::vec& Xty, const arma::vec& mu0, const arma::vec& ystar);
RcppExport SEXP _bayesDP_sigma2marginal(SEXP nSEXP, SEXP gridSEXP, SEXP XtXSEXP, SEXP SigmaBetaInvSEXP, SEXP XstarSEXP, SEXP XtySEXP, SEXP mu0SEXP, SEXP ystarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type XtX(XtXSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type SigmaBetaInv(SigmaBetaInvSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Xstar(XstarSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Xty(XtySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ystar(ystarSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma2marginal(n, grid, XtX, SigmaBetaInv, Xstar, Xty, mu0, ystar));
    return rcpp_result_gen;
END_RCPP
}
// sigma2marginalmc
SEXP sigma2marginalmc(int n, const arma::vec& grid, const arma::mat& XtX, const arma::mat& SigmaBetaInv, const arma::mat& Xstar, const arma::vec& Xty, const arma::vec& mu0, const arma::vec& ystar);
RcppExport SEXP _bayesDP_sigma2marginalmc(SEXP nSEXP, SEXP gridSEXP, SEXP XtXSEXP, SEXP SigmaBetaInvSEXP, SEXP XstarSEXP, SEXP XtySEXP, SEXP mu0SEXP, SEXP ystarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type XtX(XtXSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type SigmaBetaInv(SigmaBetaInvSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Xstar(XstarSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Xty(XtySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ystar(ystarSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma2marginalmc(n, grid, XtX, SigmaBetaInv, Xstar, Xty, mu0, ystar));
    return rcpp_result_gen;
END_RCPP
}
// betaRegSampler
arma::mat betaRegSampler(const arma::vec& sigma2, const arma::mat& XtX, const arma::mat& SigmaBetaInv, const arma::vec& mu0, const arma::vec& Xty, int nsamples);
RcppExport SEXP _bayesDP_betaRegSampler(SEXP sigma2SEXP, SEXP XtXSEXP, SEXP SigmaBetaInvSEXP, SEXP mu0SEXP, SEXP XtySEXP, SEXP nsamplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type XtX(XtXSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type SigmaBetaInv(SigmaBetaInvSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Xty(XtySEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    rcpp_result_gen = Rcpp::wrap(betaRegSampler(sigma2, XtX, SigmaBetaInv, mu0, Xty, nsamples));
    return rcpp_result_gen;
END_RCPP
}
// betaRegSamplermc
arma::mat betaRegSamplermc(const arma::vec& sigma2, const arma::mat& XtX, const arma::mat& SigmaBetaInv, const arma::vec& SigmaBetaInvID, const arma::vec& mu0, const arma::vec& Xty, int nsamples);
RcppExport SEXP _bayesDP_betaRegSamplermc(SEXP sigma2SEXP, SEXP XtXSEXP, SEXP SigmaBetaInvSEXP, SEXP SigmaBetaInvIDSEXP, SEXP mu0SEXP, SEXP XtySEXP, SEXP nsamplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type XtX(XtXSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type SigmaBetaInv(SigmaBetaInvSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type SigmaBetaInvID(SigmaBetaInvIDSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Xty(XtySEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    rcpp_result_gen = Rcpp::wrap(betaRegSamplermc(sigma2, XtX, SigmaBetaInv, SigmaBetaInvID, mu0, Xty, nsamples));
    return rcpp_result_gen;
END_RCPP
}
// ppexpV
double ppexpV(double q, const arma::vec& x, const arma::vec& cuts);
RcppExport SEXP _bayesDP_ppexpV(SEXP qSEXP, SEXP xSEXP, SEXP cutsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cuts(cutsSEXP);
    rcpp_result_gen = Rcpp::wrap(ppexpV(q, x, cuts));
    return rcpp_result_gen;
END_RCPP
}
// ppexpM
arma::vec ppexpM(double q, const arma::mat& x, const arma::vec& cuts);
RcppExport SEXP _bayesDP_ppexpM(SEXP qSEXP, SEXP xSEXP, SEXP cutsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cuts(cutsSEXP);
    rcpp_result_gen = Rcpp::wrap(ppexpM(q, x, cuts));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bayesDP_sigma2marginal", (DL_FUNC) &_bayesDP_sigma2marginal, 8},
    {"_bayesDP_sigma2marginalmc", (DL_FUNC) &_bayesDP_sigma2marginalmc, 8},
    {"_bayesDP_betaRegSampler", (DL_FUNC) &_bayesDP_betaRegSampler, 6},
    {"_bayesDP_betaRegSamplermc", (DL_FUNC) &_bayesDP_betaRegSamplermc, 7},
    {"_bayesDP_ppexpV", (DL_FUNC) &_bayesDP_ppexpV, 3},
    {"_bayesDP_ppexpM", (DL_FUNC) &_bayesDP_ppexpM, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_bayesDP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
