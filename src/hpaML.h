#ifndef hpa_hpaML_H
#define hpa_hpaML_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;

Rcpp::List hpaML(Rcpp::NumericMatrix x,
	Rcpp::NumericVector pol_degrees,
	Rcpp::NumericVector tr_left,
	Rcpp::NumericVector tr_right,
	Rcpp::LogicalVector given_ind,
	Rcpp::LogicalVector omit_ind,
	NumericVector x0,
	String cov_type,
	int bootstrap,
	bool is_parallel);

List hpaLnLOptim_List(Rcpp::NumericVector x0,
	Rcpp::NumericMatrix x_data,
	Rcpp::NumericVector pol_coefficients_ind,
	Rcpp::NumericVector pol_degrees,
	Rcpp::LogicalVector given_ind,
	Rcpp::LogicalVector omit_ind,
	Rcpp::NumericVector mean_ind,
	Rcpp::NumericVector sd_ind,
	Rcpp::NumericMatrix tr_left,
	Rcpp::NumericMatrix tr_right,
	bool is_parallel);

double hpaLnLOptim(Rcpp::NumericVector x0,
	Rcpp::NumericMatrix x_data,
	Rcpp::NumericVector pol_coefficients_ind,
	Rcpp::NumericVector pol_degrees,
	Rcpp::LogicalVector given_ind,
	Rcpp::LogicalVector omit_ind,
	Rcpp::NumericVector mean_ind,
	Rcpp::NumericVector sd_ind,
	Rcpp::NumericMatrix tr_left,
	Rcpp::NumericMatrix tr_right,
	bool is_parallel);

NumericVector hpaLnLOptim_ind(Rcpp::NumericVector x0,
	Rcpp::NumericMatrix x_data,
	Rcpp::NumericVector pol_coefficients_ind,
	Rcpp::NumericVector pol_degrees,
	Rcpp::LogicalVector given_ind,
	Rcpp::LogicalVector omit_ind,
	Rcpp::NumericVector mean_ind,
	Rcpp::NumericVector sd_ind,
	Rcpp::NumericMatrix tr_left,
	Rcpp::NumericMatrix tr_right,
	bool is_parallel);

List hpaLnLOptim_grad_List(Rcpp::NumericVector x0,
	Rcpp::NumericMatrix x_data,
	Rcpp::NumericVector pol_coefficients_ind,
	Rcpp::NumericVector pol_degrees,
	Rcpp::LogicalVector given_ind,
	Rcpp::LogicalVector omit_ind,
	Rcpp::NumericVector mean_ind,
	Rcpp::NumericVector sd_ind,
	Rcpp::NumericMatrix tr_left,
	Rcpp::NumericMatrix tr_right,
	bool is_parallel);

NumericVector hpaLnLOptim_grad(Rcpp::NumericVector x0,
	Rcpp::NumericMatrix x_data,
	Rcpp::NumericVector pol_coefficients_ind,
	Rcpp::NumericVector pol_degrees,
	Rcpp::LogicalVector given_ind,
	Rcpp::LogicalVector omit_ind,
	Rcpp::NumericVector mean_ind,
	Rcpp::NumericVector sd_ind,
	Rcpp::NumericMatrix tr_left,
	Rcpp::NumericMatrix tr_right,
	bool is_parallel);

NumericMatrix hpaLnLOptim_grad_ind(Rcpp::NumericVector x0,
	Rcpp::NumericMatrix x_data,
	Rcpp::NumericVector pol_coefficients_ind,
	Rcpp::NumericVector pol_degrees,
	Rcpp::LogicalVector given_ind,
	Rcpp::LogicalVector omit_ind,
	Rcpp::NumericVector mean_ind,
	Rcpp::NumericVector sd_ind,
	Rcpp::NumericMatrix tr_left,
	Rcpp::NumericMatrix tr_right,
	bool is_parallel);

Rcpp::NumericVector predict_hpaML(List object,
	Rcpp::NumericMatrix newdata);

Rcpp::NumericVector mecdf(NumericMatrix x);

void print_summary_hpaML(Rcpp::List x);

Rcpp::List summary_hpaML(Rcpp::List model);

double AIC_hpaML(Rcpp::List model, double k);

double logLik_hpaML(Rcpp::List model);

Rcpp::StringVector starVector(Rcpp::NumericVector p_values);

#endif
