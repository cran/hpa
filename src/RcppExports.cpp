// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// dnorm_parallel
NumericVector dnorm_parallel(NumericVector x, double mean, double sd, bool is_parallel);
RcppExport SEXP _hpa_dnorm_parallel(SEXP xSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP is_parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< bool >::type is_parallel(is_parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(dnorm_parallel(x, mean, sd, is_parallel));
    return rcpp_result_gen;
END_RCPP
}
// pnorm_parallel
NumericVector pnorm_parallel(NumericVector x, double mean, double sd, bool is_parallel);
RcppExport SEXP _hpa_pnorm_parallel(SEXP xSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP is_parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< bool >::type is_parallel(is_parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(pnorm_parallel(x, mean, sd, is_parallel));
    return rcpp_result_gen;
END_RCPP
}
// hpaBinary
List hpaBinary(Rcpp::Formula formula, DataFrame data, int K, double z_mean_fixed, double z_sd_fixed, double z_constant_fixed, bool is_z_coef_first_fixed, bool is_x0_probit, bool is_sequence, NumericVector x0, String cov_type, int boot_iter, bool is_parallel);
RcppExport SEXP _hpa_hpaBinary(SEXP formulaSEXP, SEXP dataSEXP, SEXP KSEXP, SEXP z_mean_fixedSEXP, SEXP z_sd_fixedSEXP, SEXP z_constant_fixedSEXP, SEXP is_z_coef_first_fixedSEXP, SEXP is_x0_probitSEXP, SEXP is_sequenceSEXP, SEXP x0SEXP, SEXP cov_typeSEXP, SEXP boot_iterSEXP, SEXP is_parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Formula >::type formula(formulaSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type z_mean_fixed(z_mean_fixedSEXP);
    Rcpp::traits::input_parameter< double >::type z_sd_fixed(z_sd_fixedSEXP);
    Rcpp::traits::input_parameter< double >::type z_constant_fixed(z_constant_fixedSEXP);
    Rcpp::traits::input_parameter< bool >::type is_z_coef_first_fixed(is_z_coef_first_fixedSEXP);
    Rcpp::traits::input_parameter< bool >::type is_x0_probit(is_x0_probitSEXP);
    Rcpp::traits::input_parameter< bool >::type is_sequence(is_sequenceSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< String >::type cov_type(cov_typeSEXP);
    Rcpp::traits::input_parameter< int >::type boot_iter(boot_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type is_parallel(is_parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(hpaBinary(formula, data, K, z_mean_fixed, z_sd_fixed, z_constant_fixed, is_z_coef_first_fixed, is_x0_probit, is_sequence, x0, cov_type, boot_iter, is_parallel));
    return rcpp_result_gen;
END_RCPP
}
// predict_hpaBinary
NumericVector predict_hpaBinary(List object, DataFrame newdata, bool is_prob);
RcppExport SEXP _hpa_predict_hpaBinary(SEXP objectSEXP, SEXP newdataSEXP, SEXP is_probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type object(objectSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type newdata(newdataSEXP);
    Rcpp::traits::input_parameter< bool >::type is_prob(is_probSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_hpaBinary(object, newdata, is_prob));
    return rcpp_result_gen;
END_RCPP
}
// summary_hpaBinary
List summary_hpaBinary(List object);
RcppExport SEXP _hpa_summary_hpaBinary(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(summary_hpaBinary(object));
    return rcpp_result_gen;
END_RCPP
}
// print_summary_hpaBinary
void print_summary_hpaBinary(List x);
RcppExport SEXP _hpa_print_summary_hpaBinary(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    print_summary_hpaBinary(x);
    return R_NilValue;
END_RCPP
}
// plot_hpaBinary
void plot_hpaBinary(List x);
RcppExport SEXP _hpa_plot_hpaBinary(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    plot_hpaBinary(x);
    return R_NilValue;
END_RCPP
}
// AIC_hpaBinary
double AIC_hpaBinary(List object, double k);
RcppExport SEXP _hpa_AIC_hpaBinary(SEXP objectSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type object(objectSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(AIC_hpaBinary(object, k));
    return rcpp_result_gen;
END_RCPP
}
// logLik_hpaBinary
double logLik_hpaBinary(List object);
RcppExport SEXP _hpa_logLik_hpaBinary(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(logLik_hpaBinary(object));
    return rcpp_result_gen;
END_RCPP
}
// hpaML
List hpaML(NumericMatrix x, NumericVector pol_degrees, NumericVector tr_left, NumericVector tr_right, LogicalVector given_ind, LogicalVector omit_ind, NumericVector x0, String cov_type, int boot_iter, bool is_parallel);
RcppExport SEXP _hpa_hpaML(SEXP xSEXP, SEXP pol_degreesSEXP, SEXP tr_leftSEXP, SEXP tr_rightSEXP, SEXP given_indSEXP, SEXP omit_indSEXP, SEXP x0SEXP, SEXP cov_typeSEXP, SEXP boot_iterSEXP, SEXP is_parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_degrees(pol_degreesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tr_left(tr_leftSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tr_right(tr_rightSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type given_ind(given_indSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type omit_ind(omit_indSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< String >::type cov_type(cov_typeSEXP);
    Rcpp::traits::input_parameter< int >::type boot_iter(boot_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type is_parallel(is_parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(hpaML(x, pol_degrees, tr_left, tr_right, given_ind, omit_ind, x0, cov_type, boot_iter, is_parallel));
    return rcpp_result_gen;
END_RCPP
}
// predict_hpaML
NumericVector predict_hpaML(List object, NumericMatrix newdata);
RcppExport SEXP _hpa_predict_hpaML(SEXP objectSEXP, SEXP newdataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type object(objectSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type newdata(newdataSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_hpaML(object, newdata));
    return rcpp_result_gen;
END_RCPP
}
// summary_hpaML
List summary_hpaML(List object);
RcppExport SEXP _hpa_summary_hpaML(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(summary_hpaML(object));
    return rcpp_result_gen;
END_RCPP
}
// print_summary_hpaML
void print_summary_hpaML(List x);
RcppExport SEXP _hpa_print_summary_hpaML(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    print_summary_hpaML(x);
    return R_NilValue;
END_RCPP
}
// AIC_hpaML
double AIC_hpaML(List object, double k);
RcppExport SEXP _hpa_AIC_hpaML(SEXP objectSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type object(objectSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(AIC_hpaML(object, k));
    return rcpp_result_gen;
END_RCPP
}
// logLik_hpaML
double logLik_hpaML(List object);
RcppExport SEXP _hpa_logLik_hpaML(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(logLik_hpaML(object));
    return rcpp_result_gen;
END_RCPP
}
// mecdf
NumericVector mecdf(NumericMatrix x);
RcppExport SEXP _hpa_mecdf(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(mecdf(x));
    return rcpp_result_gen;
END_RCPP
}
// dhpa
NumericVector dhpa(NumericMatrix x, NumericVector pol_coefficients, NumericVector pol_degrees, LogicalVector given_ind, LogicalVector omit_ind, NumericVector mean, NumericVector sd, bool is_parallel);
RcppExport SEXP _hpa_dhpa(SEXP xSEXP, SEXP pol_coefficientsSEXP, SEXP pol_degreesSEXP, SEXP given_indSEXP, SEXP omit_indSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP is_parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_coefficients(pol_coefficientsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_degrees(pol_degreesSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type given_ind(given_indSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type omit_ind(omit_indSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< bool >::type is_parallel(is_parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(dhpa(x, pol_coefficients, pol_degrees, given_ind, omit_ind, mean, sd, is_parallel));
    return rcpp_result_gen;
END_RCPP
}
// phpa
NumericVector phpa(NumericMatrix x, NumericVector pol_coefficients, NumericVector pol_degrees, LogicalVector given_ind, LogicalVector omit_ind, NumericVector mean, NumericVector sd, bool is_parallel);
RcppExport SEXP _hpa_phpa(SEXP xSEXP, SEXP pol_coefficientsSEXP, SEXP pol_degreesSEXP, SEXP given_indSEXP, SEXP omit_indSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP is_parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_coefficients(pol_coefficientsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_degrees(pol_degreesSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type given_ind(given_indSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type omit_ind(omit_indSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< bool >::type is_parallel(is_parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(phpa(x, pol_coefficients, pol_degrees, given_ind, omit_ind, mean, sd, is_parallel));
    return rcpp_result_gen;
END_RCPP
}
// ihpa
NumericVector ihpa(NumericMatrix x_lower, NumericMatrix x_upper, NumericVector pol_coefficients, NumericVector pol_degrees, LogicalVector given_ind, LogicalVector omit_ind, NumericVector mean, NumericVector sd, bool is_parallel);
RcppExport SEXP _hpa_ihpa(SEXP x_lowerSEXP, SEXP x_upperSEXP, SEXP pol_coefficientsSEXP, SEXP pol_degreesSEXP, SEXP given_indSEXP, SEXP omit_indSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP is_parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x_lower(x_lowerSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x_upper(x_upperSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_coefficients(pol_coefficientsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_degrees(pol_degreesSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type given_ind(given_indSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type omit_ind(omit_indSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< bool >::type is_parallel(is_parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(ihpa(x_lower, x_upper, pol_coefficients, pol_degrees, given_ind, omit_ind, mean, sd, is_parallel));
    return rcpp_result_gen;
END_RCPP
}
// ehpa
NumericVector ehpa(NumericMatrix x, NumericVector pol_coefficients, NumericVector pol_degrees, LogicalVector given_ind, LogicalVector omit_ind, NumericVector mean, NumericVector sd, NumericVector expectation_powers, bool is_parallel);
RcppExport SEXP _hpa_ehpa(SEXP xSEXP, SEXP pol_coefficientsSEXP, SEXP pol_degreesSEXP, SEXP given_indSEXP, SEXP omit_indSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP expectation_powersSEXP, SEXP is_parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_coefficients(pol_coefficientsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_degrees(pol_degreesSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type given_ind(given_indSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type omit_ind(omit_indSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type expectation_powers(expectation_powersSEXP);
    Rcpp::traits::input_parameter< bool >::type is_parallel(is_parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(ehpa(x, pol_coefficients, pol_degrees, given_ind, omit_ind, mean, sd, expectation_powers, is_parallel));
    return rcpp_result_gen;
END_RCPP
}
// etrhpa
NumericVector etrhpa(NumericMatrix tr_left, NumericMatrix tr_right, NumericVector pol_coefficients, NumericVector pol_degrees, NumericVector mean, NumericVector sd, NumericVector expectation_powers, bool is_parallel);
RcppExport SEXP _hpa_etrhpa(SEXP tr_leftSEXP, SEXP tr_rightSEXP, SEXP pol_coefficientsSEXP, SEXP pol_degreesSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP expectation_powersSEXP, SEXP is_parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type tr_left(tr_leftSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type tr_right(tr_rightSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_coefficients(pol_coefficientsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_degrees(pol_degreesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type expectation_powers(expectation_powersSEXP);
    Rcpp::traits::input_parameter< bool >::type is_parallel(is_parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(etrhpa(tr_left, tr_right, pol_coefficients, pol_degrees, mean, sd, expectation_powers, is_parallel));
    return rcpp_result_gen;
END_RCPP
}
// dtrhpa
NumericVector dtrhpa(NumericMatrix x, NumericMatrix tr_left, NumericMatrix tr_right, NumericVector pol_coefficients, NumericVector pol_degrees, LogicalVector given_ind, LogicalVector omit_ind, NumericVector mean, NumericVector sd, bool is_parallel);
RcppExport SEXP _hpa_dtrhpa(SEXP xSEXP, SEXP tr_leftSEXP, SEXP tr_rightSEXP, SEXP pol_coefficientsSEXP, SEXP pol_degreesSEXP, SEXP given_indSEXP, SEXP omit_indSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP is_parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type tr_left(tr_leftSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type tr_right(tr_rightSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_coefficients(pol_coefficientsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_degrees(pol_degreesSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type given_ind(given_indSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type omit_ind(omit_indSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< bool >::type is_parallel(is_parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(dtrhpa(x, tr_left, tr_right, pol_coefficients, pol_degrees, given_ind, omit_ind, mean, sd, is_parallel));
    return rcpp_result_gen;
END_RCPP
}
// itrhpa
NumericVector itrhpa(NumericMatrix x_lower, NumericMatrix x_upper, NumericMatrix tr_left, NumericMatrix tr_right, NumericVector pol_coefficients, NumericVector pol_degrees, LogicalVector given_ind, LogicalVector omit_ind, NumericVector mean, NumericVector sd, bool is_parallel);
RcppExport SEXP _hpa_itrhpa(SEXP x_lowerSEXP, SEXP x_upperSEXP, SEXP tr_leftSEXP, SEXP tr_rightSEXP, SEXP pol_coefficientsSEXP, SEXP pol_degreesSEXP, SEXP given_indSEXP, SEXP omit_indSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP is_parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x_lower(x_lowerSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x_upper(x_upperSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type tr_left(tr_leftSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type tr_right(tr_rightSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_coefficients(pol_coefficientsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_degrees(pol_degreesSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type given_ind(given_indSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type omit_ind(omit_indSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< bool >::type is_parallel(is_parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(itrhpa(x_lower, x_upper, tr_left, tr_right, pol_coefficients, pol_degrees, given_ind, omit_ind, mean, sd, is_parallel));
    return rcpp_result_gen;
END_RCPP
}
// dhpaDiff
NumericMatrix dhpaDiff(NumericMatrix x, NumericVector pol_coefficients, NumericVector pol_degrees, LogicalVector given_ind, LogicalVector omit_ind, NumericVector mean, NumericVector sd, String type, bool is_parallel);
RcppExport SEXP _hpa_dhpaDiff(SEXP xSEXP, SEXP pol_coefficientsSEXP, SEXP pol_degreesSEXP, SEXP given_indSEXP, SEXP omit_indSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP typeSEXP, SEXP is_parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_coefficients(pol_coefficientsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_degrees(pol_degreesSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type given_ind(given_indSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type omit_ind(omit_indSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type is_parallel(is_parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(dhpaDiff(x, pol_coefficients, pol_degrees, given_ind, omit_ind, mean, sd, type, is_parallel));
    return rcpp_result_gen;
END_RCPP
}
// ihpaDiff
NumericMatrix ihpaDiff(NumericMatrix x_lower, NumericMatrix x_upper, NumericVector pol_coefficients, NumericVector pol_degrees, LogicalVector given_ind, LogicalVector omit_ind, NumericVector mean, NumericVector sd, String type, bool is_parallel);
RcppExport SEXP _hpa_ihpaDiff(SEXP x_lowerSEXP, SEXP x_upperSEXP, SEXP pol_coefficientsSEXP, SEXP pol_degreesSEXP, SEXP given_indSEXP, SEXP omit_indSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP typeSEXP, SEXP is_parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x_lower(x_lowerSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x_upper(x_upperSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_coefficients(pol_coefficientsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_degrees(pol_degreesSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type given_ind(given_indSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type omit_ind(omit_indSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type is_parallel(is_parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(ihpaDiff(x_lower, x_upper, pol_coefficients, pol_degrees, given_ind, omit_ind, mean, sd, type, is_parallel));
    return rcpp_result_gen;
END_RCPP
}
// hpaSelection
Rcpp::List hpaSelection(Rcpp::Formula selection, Rcpp::Formula outcome, DataFrame data, int z_K, int y_K, int pol_elements, bool is_Newey, NumericVector x0, bool is_Newey_loocv, double z_sd_fixed, String cov_type, int boot_iter, bool is_parallel);
RcppExport SEXP _hpa_hpaSelection(SEXP selectionSEXP, SEXP outcomeSEXP, SEXP dataSEXP, SEXP z_KSEXP, SEXP y_KSEXP, SEXP pol_elementsSEXP, SEXP is_NeweySEXP, SEXP x0SEXP, SEXP is_Newey_loocvSEXP, SEXP z_sd_fixedSEXP, SEXP cov_typeSEXP, SEXP boot_iterSEXP, SEXP is_parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Formula >::type selection(selectionSEXP);
    Rcpp::traits::input_parameter< Rcpp::Formula >::type outcome(outcomeSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type z_K(z_KSEXP);
    Rcpp::traits::input_parameter< int >::type y_K(y_KSEXP);
    Rcpp::traits::input_parameter< int >::type pol_elements(pol_elementsSEXP);
    Rcpp::traits::input_parameter< bool >::type is_Newey(is_NeweySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< bool >::type is_Newey_loocv(is_Newey_loocvSEXP);
    Rcpp::traits::input_parameter< double >::type z_sd_fixed(z_sd_fixedSEXP);
    Rcpp::traits::input_parameter< String >::type cov_type(cov_typeSEXP);
    Rcpp::traits::input_parameter< int >::type boot_iter(boot_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type is_parallel(is_parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(hpaSelection(selection, outcome, data, z_K, y_K, pol_elements, is_Newey, x0, is_Newey_loocv, z_sd_fixed, cov_type, boot_iter, is_parallel));
    return rcpp_result_gen;
END_RCPP
}
// predict_hpaSelection
List predict_hpaSelection(List object, DataFrame newdata, std::string method, bool is_cond, bool is_outcome);
RcppExport SEXP _hpa_predict_hpaSelection(SEXP objectSEXP, SEXP newdataSEXP, SEXP methodSEXP, SEXP is_condSEXP, SEXP is_outcomeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type object(objectSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type newdata(newdataSEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< bool >::type is_cond(is_condSEXP);
    Rcpp::traits::input_parameter< bool >::type is_outcome(is_outcomeSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_hpaSelection(object, newdata, method, is_cond, is_outcome));
    return rcpp_result_gen;
END_RCPP
}
// summary_hpaSelection
List summary_hpaSelection(List object);
RcppExport SEXP _hpa_summary_hpaSelection(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(summary_hpaSelection(object));
    return rcpp_result_gen;
END_RCPP
}
// print_summary_hpaSelection
void print_summary_hpaSelection(List x);
RcppExport SEXP _hpa_print_summary_hpaSelection(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    print_summary_hpaSelection(x);
    return R_NilValue;
END_RCPP
}
// plot_hpaSelection
List plot_hpaSelection(List x, bool is_outcome);
RcppExport SEXP _hpa_plot_hpaSelection(SEXP xSEXP, SEXP is_outcomeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type is_outcome(is_outcomeSEXP);
    rcpp_result_gen = Rcpp::wrap(plot_hpaSelection(x, is_outcome));
    return rcpp_result_gen;
END_RCPP
}
// AIC_hpaSelection
double AIC_hpaSelection(List object, double k);
RcppExport SEXP _hpa_AIC_hpaSelection(SEXP objectSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type object(objectSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(AIC_hpaSelection(object, k));
    return rcpp_result_gen;
END_RCPP
}
// logLik_hpaSelection
double logLik_hpaSelection(List object);
RcppExport SEXP _hpa_logLik_hpaSelection(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(logLik_hpaSelection(object));
    return rcpp_result_gen;
END_RCPP
}
// normalMoment
NumericVector normalMoment(int k, double mean, double sd, bool return_all_moments, bool is_validation);
RcppExport SEXP _hpa_normalMoment(SEXP kSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP return_all_momentsSEXP, SEXP is_validationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< bool >::type return_all_moments(return_all_momentsSEXP);
    Rcpp::traits::input_parameter< bool >::type is_validation(is_validationSEXP);
    rcpp_result_gen = Rcpp::wrap(normalMoment(k, mean, sd, return_all_moments, is_validation));
    return rcpp_result_gen;
END_RCPP
}
// truncatedNormalMoment
NumericMatrix truncatedNormalMoment(int k, NumericVector x_lower, NumericVector x_upper, double mean, double sd, NumericVector pdf_lower, NumericVector cdf_lower, NumericVector pdf_upper, NumericVector cdf_upper, NumericVector cdf_difference, bool return_all_moments, bool is_validation, bool is_parallel);
RcppExport SEXP _hpa_truncatedNormalMoment(SEXP kSEXP, SEXP x_lowerSEXP, SEXP x_upperSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP pdf_lowerSEXP, SEXP cdf_lowerSEXP, SEXP pdf_upperSEXP, SEXP cdf_upperSEXP, SEXP cdf_differenceSEXP, SEXP return_all_momentsSEXP, SEXP is_validationSEXP, SEXP is_parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_lower(x_lowerSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_upper(x_upperSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pdf_lower(pdf_lowerSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cdf_lower(cdf_lowerSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pdf_upper(pdf_upperSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cdf_upper(cdf_upperSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cdf_difference(cdf_differenceSEXP);
    Rcpp::traits::input_parameter< bool >::type return_all_moments(return_all_momentsSEXP);
    Rcpp::traits::input_parameter< bool >::type is_validation(is_validationSEXP);
    Rcpp::traits::input_parameter< bool >::type is_parallel(is_parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(truncatedNormalMoment(k, x_lower, x_upper, mean, sd, pdf_lower, cdf_lower, pdf_upper, cdf_upper, cdf_difference, return_all_moments, is_validation, is_parallel));
    return rcpp_result_gen;
END_RCPP
}
// polynomialIndex
NumericMatrix polynomialIndex(NumericVector pol_degrees);
RcppExport SEXP _hpa_polynomialIndex(SEXP pol_degreesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pol_degrees(pol_degreesSEXP);
    rcpp_result_gen = Rcpp::wrap(polynomialIndex(pol_degrees));
    return rcpp_result_gen;
END_RCPP
}
// printPolynomial
Rcpp::String printPolynomial(NumericVector pol_degrees, NumericVector pol_coefficients);
RcppExport SEXP _hpa_printPolynomial(SEXP pol_degreesSEXP, SEXP pol_coefficientsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pol_degrees(pol_degreesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pol_coefficients(pol_coefficientsSEXP);
    rcpp_result_gen = Rcpp::wrap(printPolynomial(pol_degrees, pol_coefficients));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hpa_dnorm_parallel", (DL_FUNC) &_hpa_dnorm_parallel, 4},
    {"_hpa_pnorm_parallel", (DL_FUNC) &_hpa_pnorm_parallel, 4},
    {"_hpa_hpaBinary", (DL_FUNC) &_hpa_hpaBinary, 13},
    {"_hpa_predict_hpaBinary", (DL_FUNC) &_hpa_predict_hpaBinary, 3},
    {"_hpa_summary_hpaBinary", (DL_FUNC) &_hpa_summary_hpaBinary, 1},
    {"_hpa_print_summary_hpaBinary", (DL_FUNC) &_hpa_print_summary_hpaBinary, 1},
    {"_hpa_plot_hpaBinary", (DL_FUNC) &_hpa_plot_hpaBinary, 1},
    {"_hpa_AIC_hpaBinary", (DL_FUNC) &_hpa_AIC_hpaBinary, 2},
    {"_hpa_logLik_hpaBinary", (DL_FUNC) &_hpa_logLik_hpaBinary, 1},
    {"_hpa_hpaML", (DL_FUNC) &_hpa_hpaML, 10},
    {"_hpa_predict_hpaML", (DL_FUNC) &_hpa_predict_hpaML, 2},
    {"_hpa_summary_hpaML", (DL_FUNC) &_hpa_summary_hpaML, 1},
    {"_hpa_print_summary_hpaML", (DL_FUNC) &_hpa_print_summary_hpaML, 1},
    {"_hpa_AIC_hpaML", (DL_FUNC) &_hpa_AIC_hpaML, 2},
    {"_hpa_logLik_hpaML", (DL_FUNC) &_hpa_logLik_hpaML, 1},
    {"_hpa_mecdf", (DL_FUNC) &_hpa_mecdf, 1},
    {"_hpa_dhpa", (DL_FUNC) &_hpa_dhpa, 8},
    {"_hpa_phpa", (DL_FUNC) &_hpa_phpa, 8},
    {"_hpa_ihpa", (DL_FUNC) &_hpa_ihpa, 9},
    {"_hpa_ehpa", (DL_FUNC) &_hpa_ehpa, 9},
    {"_hpa_etrhpa", (DL_FUNC) &_hpa_etrhpa, 8},
    {"_hpa_dtrhpa", (DL_FUNC) &_hpa_dtrhpa, 10},
    {"_hpa_itrhpa", (DL_FUNC) &_hpa_itrhpa, 11},
    {"_hpa_dhpaDiff", (DL_FUNC) &_hpa_dhpaDiff, 9},
    {"_hpa_ihpaDiff", (DL_FUNC) &_hpa_ihpaDiff, 10},
    {"_hpa_hpaSelection", (DL_FUNC) &_hpa_hpaSelection, 13},
    {"_hpa_predict_hpaSelection", (DL_FUNC) &_hpa_predict_hpaSelection, 5},
    {"_hpa_summary_hpaSelection", (DL_FUNC) &_hpa_summary_hpaSelection, 1},
    {"_hpa_print_summary_hpaSelection", (DL_FUNC) &_hpa_print_summary_hpaSelection, 1},
    {"_hpa_plot_hpaSelection", (DL_FUNC) &_hpa_plot_hpaSelection, 2},
    {"_hpa_AIC_hpaSelection", (DL_FUNC) &_hpa_AIC_hpaSelection, 2},
    {"_hpa_logLik_hpaSelection", (DL_FUNC) &_hpa_logLik_hpaSelection, 1},
    {"_hpa_normalMoment", (DL_FUNC) &_hpa_normalMoment, 5},
    {"_hpa_truncatedNormalMoment", (DL_FUNC) &_hpa_truncatedNormalMoment, 13},
    {"_hpa_polynomialIndex", (DL_FUNC) &_hpa_polynomialIndex, 1},
    {"_hpa_printPolynomial", (DL_FUNC) &_hpa_printPolynomial, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_hpa(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
