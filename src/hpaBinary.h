#ifndef hpa_hpaBinary_H
#define hpa_hpaBinary_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;

Rcpp::List hpaBinary(Rcpp::Formula formula,
	Rcpp::DataFrame data,
	int K,
	double z_mean_fixed,
	double z_sd_fixed,
	double z_constant_fixed,
	bool is_z_coef_first_fixed,
	bool is_x0_probit,
	bool is_sequence,
	Rcpp::NumericVector x0,
	Rcpp::String cov_type,
	int boot_iter,
	bool is_parallel);

List hpaBinaryLnLOptim_List(Rcpp::NumericVector x0,
	      Rcpp::List is_List,
	      arma::vec z_1,
	      arma::vec z_0,
	      arma::mat z_d_1,
	      arma::mat z_d_0,
	      int K,
	      double z_mean_fixed,
	      double z_sd_fixed,
	      double z_constant_fixed0,
	      Rcpp::NumericVector pol_coefficients_ind,
	      int z_mean_ind,
	      int z_sd_ind,
	      Rcpp::NumericVector z_coef_ind,
	      bool is_parallel);

double hpaBinaryLnLOptim(Rcpp::NumericVector x0,
	      Rcpp::List is_List,
	      arma::vec z_1,
	      arma::vec z_0,
	      arma::mat z_d_1,
	      arma::mat z_d_0,
	      int K,
	      double z_mean_fixed,
	      double z_sd_fixed,
	      double z_constant_fixed0,
	      Rcpp::NumericVector pol_coefficients_ind,
	      int z_mean_ind,
	      int z_sd_ind,
	      Rcpp::NumericVector z_coef_ind,
	      bool is_parallel);

NumericVector hpaBinaryLnLOptim_ind(Rcpp::NumericVector x0,
	      Rcpp::List is_List,
	      arma::vec z_1,
	      arma::vec z_0,
	      arma::mat z_d_1,
	      arma::mat z_d_0,
	      int K,
	      double z_mean_fixed,
	      double z_sd_fixed,
	      double z_constant_fixed0,
	      Rcpp::NumericVector pol_coefficients_ind,
	      int z_mean_ind,
	      int z_sd_ind,
	      Rcpp::NumericVector z_coef_ind,
	      bool is_parallel);

Rcpp::List hpaBinaryLnLOptim_grad_List(Rcpp::NumericVector x0,
	      Rcpp::List is_List,
	      arma::vec z_1,
	      arma::vec z_0,
	      arma::mat z_d_1,
	      arma::mat z_d_0,
	      int K,
	      double z_mean_fixed,
	      double z_sd_fixed,
	      double z_constant_fixed0,
	      Rcpp::NumericVector pol_coefficients_ind,
	      int z_mean_ind,
	      int z_sd_ind,
	      Rcpp::NumericVector z_coef_ind,
	      bool is_parallel);

Rcpp::NumericVector hpaBinaryLnLOptim_grad(Rcpp::NumericVector x0,
	      Rcpp::List is_List,
	      arma::vec z_1,
	      arma::vec z_0,
	      arma::mat z_d_1,
	      arma::mat z_d_0,
	      int K,
	      double z_mean_fixed,
	      double z_sd_fixed,
	      double z_constant_fixed0,
	      Rcpp::NumericVector pol_coefficients_ind,
	      int z_mean_ind,
	      int z_sd_ind,
	      Rcpp::NumericVector z_coef_ind,
	      bool is_parallel);

Rcpp::NumericMatrix hpaBinaryLnLOptim_grad_ind(Rcpp::NumericVector x0,
	      Rcpp::List is_List,
	      arma::vec z_1,
	      arma::vec z_0,
	      arma::mat z_d_1,
	      arma::mat z_d_0,
	      int K,
	      double z_mean_fixed,
	      double z_sd_fixed,
	      double z_constant_fixed0,
	      Rcpp::NumericVector pol_coefficients_ind,
	      int z_mean_ind,
	      int z_sd_ind,
	      Rcpp::NumericVector z_coef_ind,
	      bool is_parallel);

Rcpp::NumericVector predict_hpaBinary(Rcpp::List object, 
	Rcpp::DataFrame newdata, 
	bool is_prob);

Rcpp::List summary_hpaBinary(Rcpp::List object);

void print_summary_hpaBinary(Rcpp::List x);

void plot_hpaBinary(Rcpp::List x);

double AIC_hpaBinary(Rcpp::List object, double k);

double logLik_hpaBinary(Rcpp::List object);

#endif
