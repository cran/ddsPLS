// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// modelddsPLSCpp_Rcpp
Rcpp::List modelddsPLSCpp_Rcpp(const Eigen::MatrixXd U, const Eigen::MatrixXd V, const Eigen::MatrixXd X, const Eigen::MatrixXd Y, const Eigen::VectorXd lambdas, const int R, const int n, const int p, const int q, const Eigen::VectorXd lambda0);
RcppExport SEXP _ddsPLS_modelddsPLSCpp_Rcpp(SEXP USEXP, SEXP VSEXP, SEXP XSEXP, SEXP YSEXP, SEXP lambdasSEXP, SEXP RSEXP, SEXP nSEXP, SEXP pSEXP, SEXP qSEXP, SEXP lambda0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type U(USEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type V(VSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< const int >::type R(RSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type q(qSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type lambda0(lambda0SEXP);
    rcpp_result_gen = Rcpp::wrap(modelddsPLSCpp_Rcpp(U, V, X, Y, lambdas, R, n, p, q, lambda0));
    return rcpp_result_gen;
END_RCPP
}
// bootstrap_Rcpp
Rcpp::List bootstrap_Rcpp(const Eigen::MatrixXd U, const Eigen::MatrixXd V, const Eigen::MatrixXd X, const Eigen::MatrixXd Y, const Eigen::VectorXd lambdas, const Eigen::VectorXd lambda_prev, const int R, const int n_B, const bool doBoot, const int n, const int p, const int q, const int N_lambdas, const Eigen::VectorXd lambda0);
RcppExport SEXP _ddsPLS_bootstrap_Rcpp(SEXP USEXP, SEXP VSEXP, SEXP XSEXP, SEXP YSEXP, SEXP lambdasSEXP, SEXP lambda_prevSEXP, SEXP RSEXP, SEXP n_BSEXP, SEXP doBootSEXP, SEXP nSEXP, SEXP pSEXP, SEXP qSEXP, SEXP N_lambdasSEXP, SEXP lambda0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type U(USEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type V(VSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type lambda_prev(lambda_prevSEXP);
    Rcpp::traits::input_parameter< const int >::type R(RSEXP);
    Rcpp::traits::input_parameter< const int >::type n_B(n_BSEXP);
    Rcpp::traits::input_parameter< const bool >::type doBoot(doBootSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type q(qSEXP);
    Rcpp::traits::input_parameter< const int >::type N_lambdas(N_lambdasSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type lambda0(lambda0SEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrap_Rcpp(U, V, X, Y, lambdas, lambda_prev, R, n_B, doBoot, n, p, q, N_lambdas, lambda0));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ddsPLS_modelddsPLSCpp_Rcpp", (DL_FUNC) &_ddsPLS_modelddsPLSCpp_Rcpp, 10},
    {"_ddsPLS_bootstrap_Rcpp", (DL_FUNC) &_ddsPLS_bootstrap_Rcpp, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_ddsPLS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
