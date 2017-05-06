// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cpp_calc_critical
double cpp_calc_critical(int r, NumericVector lower, NumericVector upper, double error_spend, NumericVector information, double theta, String side);
RcppExport SEXP gscounts_cpp_calc_critical(SEXP rSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP error_spendSEXP, SEXP informationSEXP, SEXP thetaSEXP, SEXP sideSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< double >::type error_spend(error_spendSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type information(informationSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< String >::type side(sideSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_calc_critical(r, lower, upper, error_spend, information, theta, side));
    return rcpp_result_gen;
END_RCPP
}
// cpp_pmultinorm
double cpp_pmultinorm(int r, NumericVector lower, NumericVector upper, NumericVector information, double theta);
RcppExport SEXP gscounts_cpp_pmultinorm(SEXP rSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP informationSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type information(informationSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_pmultinorm(r, lower, upper, information, theta));
    return rcpp_result_gen;
END_RCPP
}