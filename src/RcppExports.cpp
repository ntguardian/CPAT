// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// besselJ_zeros_cpp
NumericVector besselJ_zeros_cpp(const double& nu, const unsigned& a, const unsigned& b);
RcppExport SEXP _CPAT_besselJ_zeros_cpp(SEXP nuSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(besselJ_zeros_cpp(nu, a, b));
    return rcpp_result_gen;
END_RCPP
}
// stat_Vn_cpp
List stat_Vn_cpp(const NumericVector& dat, const double& kn, const double& tau, const bool& use_kernel_var, const NumericVector& lrv_est, const bool& get_all_vals);
RcppExport SEXP _CPAT_stat_Vn_cpp(SEXP datSEXP, SEXP knSEXP, SEXP tauSEXP, SEXP use_kernel_varSEXP, SEXP lrv_estSEXP, SEXP get_all_valsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type dat(datSEXP);
    Rcpp::traits::input_parameter< const double& >::type kn(knSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const bool& >::type use_kernel_var(use_kernel_varSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lrv_est(lrv_estSEXP);
    Rcpp::traits::input_parameter< const bool& >::type get_all_vals(get_all_valsSEXP);
    rcpp_result_gen = Rcpp::wrap(stat_Vn_cpp(dat, kn, tau, use_kernel_var, lrv_est, get_all_vals));
    return rcpp_result_gen;
END_RCPP
}
// stat_Zn_cpp
List stat_Zn_cpp(const NumericVector& dat, const double& kn, const bool& use_kernel_var, const NumericVector& lrv_est, const bool& get_all_vals);
RcppExport SEXP _CPAT_stat_Zn_cpp(SEXP datSEXP, SEXP knSEXP, SEXP use_kernel_varSEXP, SEXP lrv_estSEXP, SEXP get_all_valsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type dat(datSEXP);
    Rcpp::traits::input_parameter< const double& >::type kn(knSEXP);
    Rcpp::traits::input_parameter< const bool& >::type use_kernel_var(use_kernel_varSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lrv_est(lrv_estSEXP);
    Rcpp::traits::input_parameter< const bool& >::type get_all_vals(get_all_valsSEXP);
    rcpp_result_gen = Rcpp::wrap(stat_Zn_cpp(dat, kn, use_kernel_var, lrv_est, get_all_vals));
    return rcpp_result_gen;
END_RCPP
}
// stat_Zn_reg_cpp
List stat_Zn_reg_cpp(const NumericMatrix& X_input, const NumericVector& y_input, const double& kn, NumericVector lrv_est, const bool& get_all_vals, const bool& fast);
RcppExport SEXP _CPAT_stat_Zn_reg_cpp(SEXP X_inputSEXP, SEXP y_inputSEXP, SEXP knSEXP, SEXP lrv_estSEXP, SEXP get_all_valsSEXP, SEXP fastSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X_input(X_inputSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y_input(y_inputSEXP);
    Rcpp::traits::input_parameter< const double& >::type kn(knSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lrv_est(lrv_estSEXP);
    Rcpp::traits::input_parameter< const bool& >::type get_all_vals(get_all_valsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type fast(fastSEXP);
    rcpp_result_gen = Rcpp::wrap(stat_Zn_reg_cpp(X_input, y_input, kn, lrv_est, get_all_vals, fast));
    return rcpp_result_gen;
END_RCPP
}
// get_lrv_vec_cpp
NumericVector get_lrv_vec_cpp(const NumericMatrix& Y, const NumericVector& kern, const int& max_l);
RcppExport SEXP _CPAT_get_lrv_vec_cpp(SEXP YSEXP, SEXP kernSEXP, SEXP max_lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type kern(kernSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_l(max_lSEXP);
    rcpp_result_gen = Rcpp::wrap(get_lrv_vec_cpp(Y, kern, max_l));
    return rcpp_result_gen;
END_RCPP
}
// cond_var_gradient_hessian_cpp
List cond_var_gradient_hessian_cpp(const NumericVector& var, const NumericVector& eps, const double& omega, const double& alpha, const double& beta, const NumericVector& init_vals);
RcppExport SEXP _CPAT_cond_var_gradient_hessian_cpp(SEXP varSEXP, SEXP epsSEXP, SEXP omegaSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP init_valsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type var(varSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const double& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type init_vals(init_valsSEXP);
    rcpp_result_gen = Rcpp::wrap(cond_var_gradient_hessian_cpp(var, eps, omega, alpha, beta, init_vals));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CPAT_besselJ_zeros_cpp", (DL_FUNC) &_CPAT_besselJ_zeros_cpp, 3},
    {"_CPAT_stat_Vn_cpp", (DL_FUNC) &_CPAT_stat_Vn_cpp, 6},
    {"_CPAT_stat_Zn_cpp", (DL_FUNC) &_CPAT_stat_Zn_cpp, 5},
    {"_CPAT_stat_Zn_reg_cpp", (DL_FUNC) &_CPAT_stat_Zn_reg_cpp, 6},
    {"_CPAT_get_lrv_vec_cpp", (DL_FUNC) &_CPAT_get_lrv_vec_cpp, 3},
    {"_CPAT_cond_var_gradient_hessian_cpp", (DL_FUNC) &_CPAT_cond_var_gradient_hessian_cpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_CPAT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
