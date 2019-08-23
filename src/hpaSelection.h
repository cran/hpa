#ifndef hpa_hpaSelection_H
#define hpa_hpaSelection_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;

Rcpp::List hpaSelection(Rcpp::Formula selection,
	Rcpp::Formula outcome,
	Rcpp::DataFrame data,
	int z_K,
	int y_K,
	int pol_elements,
	bool is_Newey,
	Rcpp::NumericVector x0);

double hpaSelectionLnLOptim(Rcpp::NumericVector x0,
	Rcpp::List ind_List,
	arma::vec y_1,
	arma::vec y_0,
	arma::vec z_1,
	arma::vec z_0,
	arma::mat y_d_1,
	arma::mat y_d_0,
	arma::mat z_d_1,
	arma::mat z_d_0,
	Rcpp::NumericVector pol_degrees,
	double z_sd,
	bool is_minus);

Rcpp::List predict_hpaSelection(Rcpp::List object, 
	Rcpp::DataFrame newdata, 
	std::string method, 
	bool is_cond,
	bool is_outcome);

double AIC_hpaSelection(Rcpp::List model, double k);

double logLik_hpaSelection(Rcpp::List model);

Rcpp::List plot_hpaSelection(Rcpp::List x, bool is_outcome);

void print_summary_hpaSelection(Rcpp::List x);

Rcpp::List summary_hpaSelection(Rcpp::List object);

#endif
