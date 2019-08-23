#ifndef hpa_hpaMain_H
#define hpa_hpaMain_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;
	
NumericVector hpaMain(
	NumericMatrix x_lower,
	NumericMatrix x_upper,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	String type,
	LogicalVector given_ind,
	LogicalVector omit_ind,
	NumericVector mean,
	NumericVector sd,
	NumericVector expectation_powers,
	NumericMatrix pdf_lower,
	NumericMatrix cdf_lower,
	NumericMatrix pdf_upper,
	NumericMatrix cdf_upper,
	NumericMatrix cdf_difference);

NumericVector dhpa(
	NumericMatrix x,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	LogicalVector given_ind,
	LogicalVector omit_ind,
	NumericVector mean,
	NumericVector sd);

NumericVector phpa(
	NumericMatrix x,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	LogicalVector given_ind,
	LogicalVector omit_ind,
	NumericVector mean,
	NumericVector sd);

NumericVector ihpa(
	NumericMatrix x_lower,
	NumericMatrix x_upper,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	LogicalVector given_ind,
	LogicalVector omit_ind,
	NumericVector mean,
	NumericVector sd);

NumericVector ehpa(NumericMatrix x,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	LogicalVector given_ind,
	LogicalVector omit_ind,
	NumericVector mean,
	NumericVector sd,
	NumericVector expectation_powers);

//Truncated distribution

NumericVector etrhpa(
	NumericMatrix tr_left,
	NumericMatrix tr_right,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	NumericVector mean,
	NumericVector sd,
	NumericVector expectation_powers);

NumericVector dtrhpa(
	NumericMatrix x,
	NumericMatrix tr_left,
	NumericMatrix tr_right,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	LogicalVector given_ind,
	LogicalVector omit_ind,
	NumericVector mean,
	NumericVector sd);

NumericVector itrhpa(
	NumericMatrix x_lower,
	NumericMatrix x_upper,
	NumericMatrix x_left,
	NumericMatrix x_right,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	LogicalVector given_ind,
	LogicalVector omit_ind,
	NumericVector mean,
	NumericVector sd);

#endif
