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
	Rcpp::NumericVector x0,
	Rcpp::String cov_type,
	int bootstrap,
	bool is_parallel,
	Rcpp::String opt_type,
	Rcpp::List opt_control);

List hpaLnLOptim_List(Rcpp::NumericVector x0, Rcpp::List hpaML_args);

double hpaLnLOptim(Rcpp::NumericVector x0, Rcpp::List hpaML_args);

NumericVector hpaLnLOptim_ind(Rcpp::NumericVector x0, Rcpp::List hpaML_args);

List hpaLnLOptim_grad_List(Rcpp::NumericVector x0, Rcpp::List hpaML_args);

NumericVector hpaLnLOptim_grad(Rcpp::NumericVector x0, Rcpp::List hpaML_args);

NumericMatrix hpaLnLOptim_grad_ind(Rcpp::NumericVector x0, Rcpp::List hpaML_args);

NumericMatrix hpaLnLOptim_hessian(Rcpp::NumericVector x0, Rcpp::List hpaML_args);

Rcpp::NumericVector predict_hpaML(Rcpp::List object,
	Rcpp::NumericMatrix newdata);

Rcpp::NumericVector mecdf(NumericMatrix x);

void print_summary_hpaML(Rcpp::List x);

Rcpp::List summary_hpaML(Rcpp::List model);

double AIC_hpaML(Rcpp::List model, double k);

double logLik_hpaML(Rcpp::List model);

Rcpp::StringVector starVector(Rcpp::NumericVector p_values);

#endif
