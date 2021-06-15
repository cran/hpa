#ifndef hpa_hpaMain_H
#define hpa_hpaMain_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;
	
List hpaMain(
	NumericVector x_lower_vec,
	NumericVector x_upper_vec,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	String type,
	NumericVector given_ind,
	NumericVector omit_ind,
	NumericVector mean,
	NumericVector sd,
	NumericVector expectation_powers,
	NumericMatrix pdf_lower,
	NumericMatrix cdf_lower,
	NumericMatrix pdf_upper,
	NumericMatrix cdf_upper,
	NumericMatrix cdf_difference,
	String grad_type,
	bool is_parallel,
	bool log,
	bool is_validation);

NumericVector dhpa(
	NumericVector x,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	NumericVector given_ind,
	NumericVector omit_ind,
	NumericVector mean,
	NumericVector sd,
	bool is_parallel,
	bool log,
	bool is_validation);

NumericVector phpa(
	NumericVector x,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	NumericVector given_ind,
	NumericVector omit_ind,
	NumericVector mean,
	NumericVector sd,
	bool is_parallel,
	bool log,
	bool is_validation);

NumericVector ihpa(
	NumericVector x_lower,
	NumericVector x_upper,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	NumericVector given_ind,
	NumericVector omit_ind,
	NumericVector mean,
	NumericVector sd,
	bool is_parallel,
	bool log,
	bool is_validation);

NumericVector ehpa(NumericVector x,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	NumericVector given_ind,
	NumericVector omit_ind,
	NumericVector mean,
	NumericVector sd,
	NumericVector expectation_powers,
	bool is_parallel,
	bool is_validation);

//Truncated distribution

NumericVector etrhpa(
	NumericVector tr_left,
	NumericVector tr_right,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	NumericVector mean,
	NumericVector sd,
	NumericVector expectation_powers,
	bool is_parallel,
	bool is_validation);

NumericVector dtrhpa(
	NumericVector x,
	NumericVector tr_left,
	NumericVector tr_right,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	NumericVector given_ind,
	NumericVector omit_ind,
	NumericVector mean,
	NumericVector sd,
	bool is_parallel,
	bool log,
	bool is_validation);

NumericVector itrhpa(
	NumericVector x_lower,
	NumericVector x_upper,
	NumericVector x_left,
	NumericVector x_right,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	NumericVector given_ind,
	NumericVector omit_ind,
	NumericVector mean,
	NumericVector sd,
	bool is_parallel,
	bool log,
	bool is_validation);

NumericMatrix dhpaDiff(
	NumericVector x,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	NumericVector given_ind,
	NumericVector omit_ind,
	NumericVector mean,
	NumericVector sd,
	String type,
	bool is_parallel,
	bool log,
	bool is_validation);

NumericMatrix ehpaDiff(
    NumericVector x,
    NumericVector pol_coefficients,
    NumericVector pol_degrees,
    NumericVector given_ind,
    NumericVector omit_ind,
    NumericVector mean,
    NumericVector sd,
    NumericVector expectation_powers,
    String type,
    bool is_parallel,
    bool log,
    bool is_validation);

NumericMatrix dehpaDiff(
    NumericVector x,
    NumericVector pol_coefficients,
    NumericVector pol_degrees,
    NumericVector given_ind,
    NumericVector omit_ind,
    NumericVector mean,
    NumericVector sd,
    NumericVector expectation_powers,
    String diffType,
    String type,
    bool is_parallel,
    bool log,
    bool is_validation);

NumericMatrix ihpaDiff(
	NumericVector x_lower,
	NumericVector x_upper,
	NumericVector pol_coefficients,
	NumericVector pol_degrees,
	NumericVector given_ind,
	NumericVector omit_ind,
	NumericVector mean,
	NumericVector sd,
	String type,
	bool is_parallel,
	bool log,
	bool is_validation);

NumericVector qhpa(
    NumericVector p,
    NumericMatrix x,
    NumericVector pol_coefficients,
    NumericVector pol_degrees,
    NumericVector given_ind,
    NumericVector omit_ind,
    NumericVector mean,
    NumericVector sd);

double qhpa_opt(
    NumericVector par,
    NumericVector x,
    NumericVector p,
    NumericVector pol_coefficients,
    NumericVector pol_degrees,
    NumericVector given_ind,
    NumericVector omit_ind,
    NumericVector mean,
    NumericVector sd,
    int pol_degrees_n,
    int x_ind);

NumericMatrix rhpa(
    int n,
    NumericVector pol_coefficients,
    NumericVector pol_degrees,
    NumericVector mean,
    NumericVector sd);

#endif
