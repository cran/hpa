#include "hpaMain.h"
#include "hpaML.h"
#include "hpaSelection.h"
#include "hpaBinary.h"
#include "polynomialIndex.h"

#include <RcppArmadillo.h>

using namespace RcppArmadillo;

// [[Rcpp::depends(RcppArmadillo)]]


//' Perform semi-nonparametric selection model estimation
//' @description This function performs semi-nonparametric selection model estimation
//' via hermite polynomial densities approximation.
//' @param selection an object of class "formula" (or one that can be coerced to that class): 
//' a symbolic description of the selection equation form.
//' @param outcome an object of class "formula" (or one that can be coerced to that class): 
//' a symbolic description of the outcome equation form.
//' @template data_Template
//' @template z_K_Template
//' @template y_K_Template
//' @param pol_elements number of conditional expectation approximating terms for Newey method.
//' @param is_Newey logical; if TRUE then returns only Newey's method estimation results (default value is FALSE).
//' @template x0_selection_Template
//' @template hpa_likelihood_details_Template
//' @template GN_details_Template
//' @template first_coef_Template
//' @details Note that coefficient for the first
//' independent variable in \code{selection} will be fixed to 1.
//' @template sd_adjust_Template
//' @template is_numeric_selection_Template
//' @template parametric_paradigm_Template
//' @template Newey_details_Template
//' @details Note that selection equation dependent variables should have exactly two levels (0 and 1) where "0" states for the selection results 
//' which leads to unobservable values of dependent variable in outcome equation.
//' @template Mroz_reference_Template
//' @template optim_details_Template
//' @return This function returns an object of class "hpaSelection".\cr \cr
//' An object of class "hpaSelection" is a list containing the following components:
//' \itemize{
//' \item \code{optim} - \code{\link[stats]{optim}} function output.
//' \item \code{x1} - numeric vector of distribution parameters estimates.
//' \item \code{Newey} - list containing information concerning Newey's method estimation results.
//' \item \code{z_mean} - estimate of the hermite polynomial mean parameter related to selection equation random error marginal distribution.
//' \item \code{y_mean} - estimate of the hermite polynomial mean parameter related to outcome equation random error marginal distribution.
//' \item \code{z_sd} - adjusted value of sd parameter related to selection equation random error marginal distribution.
//' \item \code{y_sd} - estimate of the hermite polynomial sd parameter related to outcome equation random error marginal distribution.
//' \item \code{pol_coefficients} - polynomial coefficients estimates.
//' \item \code{pol_degrees} - numeric vector which first element is \code{z_K} and the second is \code{y_K}.
//' \item \code{z_coef} - selection equation regression coefficients estimates.
//' \item \code{y_coef} - outcome equation regression coefficients estimates.
//' \item \code{cov_matrix} - estimated parameters covariance matrix estimate.
//' \item \code{results} - numeric matrix representing estimation results.
//' \item \code{log-likelihood} - value of Log-Likelihood function.
//' \item \code{AIC} - AIC value.
//' \item \code{re_moments} - list which contains information about random errors expectations, variances and correlation.
//' \item \code{data_List} - list containing model variables and their partiotion according to outcome and selection equations.
//' \item \code{n_obs} - number of observations.
//' \item \code{ind_List} - list which contains information about parameters indexes in \code{x1}.
//' \item \code{selection_formula} - the same as \code{selection} input parameter.
//' \item \code{outcome_formula} - the same as \code{outcome} input parameter.}
//' Abovementioned list \code{Newey} has class "hpaNewey" and contains the following components:
//' \itemize{
//' \item \code{y_coef} - regression coefficients estimates (except constant term which is part of conditional expectation approximating polynomial).
//' \item \code{z_coef} - regression coefficients estimates related to selection equation.
//' \item \code{constant_biased} - biased estimate of constant term.
//' \item \code{inv_mills} - inverse mills rations estimates and their powers (including constant).
//' \item \code{inv_mills_coef} - coefficients related to \code{inv_mills}.
//' \item \code{pol_elements} - the same as \code{pol_elements} input parameter.
//' \item \code{outcome_exp_cond} - dependend variable conditional expectation estimates.
//' \item \code{selection_exp} - selection equation random error expectation estimate.
//' \item \code{selection_var} - selection equation random error variance estimate.
//' \item \code{hpaBinaryModel} - object of class "hpaBinary" which contains selection equation estimation results.}
//' Abovementioned list \code{re_moments} contains the following components:
//' \itemize{
//' \item \code{selection_exp} - selection equation random errors expectation estimate.
//' \item \code{selection_var} - selection equation random errors variance estimate.
//' \item \code{outcome_exp} - outcome equation random errors expectation estimate.
//' \item \code{outcome_var} - outcome equation random errors variance estimate.
//' \item \code{errors_covariance} - outcome and selection equation random errors covariance estimate.
//' \item \code{rho} - outcome and selection equation random errors correlation estimate.}
//' @seealso \link[hpa]{summary.hpaSelection}, \link[hpa]{predict.hpaSelection}, \link[hpa]{plot.hpaSelection}, \link[hpa]{AIC.hpaSelection}, \link[hpa]{logLik.hpaSelection}
//' @template hpaSelection_examples_Template
//' @export	
// [[Rcpp::export]]
Rcpp::List hpaSelection(Rcpp::Formula selection,
	Rcpp::Formula outcome,
	DataFrame data,
	int z_K = 1,
	int y_K = 1,
	int pol_elements = 3,
	bool is_Newey = false,
	NumericVector x0 = NumericVector(0)) {

	//Load additional environments

	//Add in future
	//Rcpp::Environment stats_env("package:optimParallel");
	//Rcpp::Function optim = stats_env["optimParallel"];

	Rcpp::Environment stats_env("package:stats");
	Rcpp::Function optim = stats_env["optim"];
	Rcpp::Function model_frame = stats_env["model.frame"];
	Rcpp::Function na_omit_R = stats_env["na.omit"];
	Rcpp::Function na_pass = stats_env["na.pass"];
	Rcpp::Function complete_cases = stats_env["complete.cases"];

	Rcpp::Environment base_env("package:base");
	Rcpp::Function solve = base_env["solve"];
	Rcpp::Function class_R = base_env["class"];
	Rcpp::Function c_R = base_env["c"];
	Rcpp::Function cbind_R = base_env["cbind"];
	Rcpp::Function subset_R = base_env["subset"];

	//Remove NA values

	//data = na_omit_R(data);

	//Initialize polynomial structure related values

	int pol_coefficients_n = (z_K + 1) * (y_K + 1) - 1; //-1 because of restriction a(0...0)=1

	NumericVector pol_degrees = {(double)z_K, (double)y_K};

	NumericMatrix polIndex_mat = polynomialIndex(pol_degrees);

	//Working with Data

		//Extract dataframe from formula

	DataFrame z_df = model_frame(Rcpp::_["formula"] = selection, Rcpp::_["data"] = data, Rcpp::_["na.action"] = na_pass);
	DataFrame y_df = model_frame(Rcpp::_["formula"] = outcome, Rcpp::_["data"] = data, Rcpp::_["na.action"] = na_pass);

	DataFrame z_y_df = cbind_R(z_df, y_df);

	LogicalVector is_z_y_df_complete = complete_cases(z_y_df);
	LogicalVector is_z_df_complete = complete_cases(z_df);
	NumericVector z_temporal = z_y_df[0];
	LogicalVector is_y_unobs = (z_temporal == 0);

	LogicalVector df_cond = is_z_y_df_complete | (is_z_df_complete & is_y_unobs);

	z_df = subset_R(Rcpp::_["x"] = z_df, Rcpp::_["subset"] = df_cond);
	y_df = subset_R(Rcpp::_["x"] = y_df, Rcpp::_["subset"] = df_cond);

	CharacterVector z_df_names = z_df.names();
	CharacterVector y_df_names = y_df.names();

	int z_df_n = z_df.size();
	int y_df_n = y_df.size();

	DataFrame y_df_no_y = y_df;

		//Extract dependend variables

	NumericVector z = z_df[0]; //it is reference
	NumericVector y = y_df[0]; //it is reference

	int n = z.size();

	//Extract independend variable

	NumericMatrix z_d(n, z_df_n - 1);//-1 because there is no constant term
	NumericMatrix y_d(n, y_df_n - 1);//-1 because there is no constant term

	int z_d_col = z_d.ncol();
	int y_d_col = y_d.ncol();

	for (int i = 0; i < z_d_col; i++)
	{
		z_d(_, i) = NumericVector(z_df[i+1]);
	}
	for (int i = 0; i < y_d_col; i++)
	{
		y_d(_, i) = NumericVector(y_df[i + 1]);
	}

	//Create initial values vector
	bool x0_given = true;

	if (x0.size() == 0)
	{
		x0_given = false;
		x0 = NumericVector(pol_coefficients_n + 2 + z_d_col + y_d_col); //+2 for mean and sd
	}

	//Assign indexes

		//Initialize additional index and upper value for some loops

	int lower_ind = 0;

	int upper_ind = pol_coefficients_n;

		//for polynomial coefficients

	NumericVector pol_coefficients_ind(pol_coefficients_n);

	for (int i = lower_ind; i < upper_ind; i++)
	{
		pol_coefficients_ind[i] = i;
	}

		//for mean vector
	
	int z_mean_ind = pol_coefficients_n;
	int y_mean_ind = z_mean_ind + 1;

		//for sd vector

	int y_sd_ind;

	y_sd_ind = y_mean_ind + 1;

		//for z coefficients

	NumericVector z_coef_ind(z_d_col - 1);

	lower_ind = y_sd_ind + 1;

	upper_ind = lower_ind + z_d_col - 1;

	for (int i = lower_ind; i < upper_ind; i++)
	{
		z_coef_ind[i - lower_ind] = i;
	}

		//for y coefficients

	NumericVector y_coef_ind(y_d_col);

	lower_ind = upper_ind;

	upper_ind = lower_ind + y_d_col;

	for (int i = lower_ind; i < upper_ind; i++)
	{
		y_coef_ind[i - lower_ind] = i;
	}

	//Convert to arma

	arma::vec y_arma = as<arma::vec>(y);

	arma::vec z_arma = as<arma::vec>(z);

	arma::mat y_d_arma = as<arma::mat>(y_d);

	arma::mat z_d_arma = as<arma::mat>(z_d);

	//Divide into observable and unobservable samples

		//observable

	arma::vec y_1 = as<arma::vec>(y[z == 1]);
	
	int n_1 = y_1.size();

	arma::mat y_d_1 = (as<arma::mat>(y_d)).rows(arma::find(z_arma == 1));

	arma::vec z_1 = as<arma::vec>(z[z == 1]);

	arma::mat z_d_1 = (as<arma::mat>(z_d)).rows(arma::find(z_arma == 1));

		//unobservable

	arma::vec y_0 = as<arma::vec>(y[z == 0]);

	arma::mat y_d_0 = (as<arma::mat>(y_d)).rows(arma::find(z_arma == 0));

	arma::vec z_0 = as<arma::vec>(z[z == 0]);

	arma::mat z_d_0 = (as<arma::mat>(z_d)).rows(arma::find(z_arma == 0));

	List Newey;
	double z_sd = 1;

	//Get initial values from hpaBinary and Newey
	if (!x0_given | is_Newey)
	{
		//for hpaBinary
		List modelBinary;
		try
		{
			modelBinary = hpaBinary(selection,
				data, z_K,
				NA_REAL, 1, 0,
				true, true, true, NumericVector(0));
		}
		catch (std::exception &ex)
		{
			warning("Can't get initial values from binary choice model");
			forward_exception_to_r(ex);
		}

		List modelBinary_K = modelBinary[z_K];

		NumericVector z_pol_coef_temporal = modelBinary_K["pol_coefficients"];

		int z_pol_ind = 1;
		for (int i = 1; i < pol_coefficients_n; i++)
		{
			if (polIndex_mat(1, i) == 0)
			{
				x0[pol_coefficients_ind[i - 1]] = z_pol_coef_temporal[z_pol_ind];
				z_pol_ind++;
			}
		}

		NumericVector z_coef_temporal = modelBinary_K["coefficients"];
		NumericVector z_coef_ind_temporal = z_coef_temporal[Rcpp::Range(1, z_d_col - 1)];
		x0[z_coef_ind] = z_coef_ind_temporal;
		x0[z_mean_ind] = modelBinary_K["mean"];
		z_sd = modelBinary_K["sd"];

		NumericVector z_latent = wrap(z_d_arma * as<arma::vec>(z_coef_temporal));
		double z_exp = modelBinary_K["errors_exp"];
		double z_var = modelBinary_K["errors_var"];
		z_latent = (z_latent - z_exp) / sqrt(z_var); //standartize for mills ratio

		NumericVector z_latent_1 = z_latent[z == 1];
		NumericVector z_latent_0 = z_latent[z == 0];

		//Newey with 3 approximating polynomial elements

			//prepare data

		arma::mat y_d_1_Newey = arma::mat(y_1.size(), y_d_col + pol_elements + 1, arma::fill::ones);
		
		NumericMatrix z_mills = NumericMatrix(n_1, pol_elements + 1);

		for (int i = 0; i < y_d_col; i++)
		{
			y_d_1_Newey.col(i) = y_d_1.col(i);
		}
		for (int i = y_d_col; i < (y_d_col + pol_elements + 1); i++)
		{
			NumericVector z_mills_temporal = dnorm(z_latent_1) / pnorm(z_latent_1);
			z_mills(_, i-y_d_col) = pow(z_mills_temporal, i - y_d_col);
			z_mills_temporal = z_mills(_, i - y_d_col);
			y_d_1_Newey.col(i) = as<arma::vec>(z_mills_temporal);
		}

		//estimate coefficients

		arma::mat y_coef_Newey_arma;
		y_coef_Newey_arma = inv(y_d_1_Newey.t() * y_d_1_Newey) * y_d_1_Newey.t() * y_1;

		NumericVector y_coef_Newey = wrap(y_coef_Newey_arma);

		//assign coefficients to x0

		NumericVector y_coef_Newey_temporal = y_coef_Newey[Rcpp::Range(0, y_d_col - 1)];

		x0[y_coef_ind] = y_coef_Newey_temporal;

		x0[y_mean_ind] = y_coef_Newey[y_d_col];

		arma::mat predictions_ls = y_d_1_Newey * y_coef_Newey_arma - y_1;
		arma::mat residuals = predictions_ls.t() * predictions_ls;
		x0[y_sd_ind] = sqrt(residuals(0, 0) / (n_1 - 1 - y_d_col));

		//get additional values for inverse mills ratios

		NumericVector y_coef_Newey_mills = y_coef_Newey[Rcpp::Range(y_d_col, y_d_col + pol_elements)];

		NumericVector z_expect = NumericVector(n_1);

		for (int i = 0; i < (pol_elements + 1); i++)
		{
			z_expect = z_expect + y_coef_Newey_mills[i] * z_mills(_, i);
		}

		//summarize output for Newey

		Newey = List::create(Named("y_coef") = y_coef_Newey_temporal,
			Named("z_coef") = z_coef_ind_temporal,
			Named("constant_biased") = x0[y_mean_ind], 
			Named("sd_biased") = x0[y_sd_ind],
			Named("inv_mills") = z_mills,
			Named("inv_mills_coef") = y_coef_Newey_mills,
			Named("pol_elements") = pol_elements,
			Named("outcome_exp_cond") = z_expect,
			Named("selection_exp") = z_exp,
			Named("selection_var") = z_var,
			Named("z_mean") = x0[z_mean_ind],
			Named("z_sd") = z_sd,
			Named("y_mean") = x0[y_mean_ind],
			Named("y_sd") = x0[y_sd_ind],
			Named("hpaBinaryModel") = modelBinary);

		Newey.attr("class") = "hpaNewey";

		if (is_Newey)
		{
			return(Newey);
		}
	}

	//Create list for some variables because unfortunatelly optim function has limitation for the
	//parameters number (parameters itself not estimated)

	List ind_List = List::create(Named("pol_coefficients_ind") = pol_coefficients_ind,
		Named("z_mean_ind") = z_mean_ind,
		Named("y_mean_ind") = y_mean_ind,
		Named("y_sd_ind") = y_sd_ind,
		Named("y_coef_ind") = y_coef_ind,
		Named("z_coef_ind") = z_coef_ind
	);

	//Apply optimization routine

	List PGN_control = List::create(Named("maxit") = 1000000);

	List optim_results = optim(
		Rcpp::_["par"] = x0,
		Rcpp::_["fn"] = Rcpp::InternalFunction(&hpaSelectionLnLOptim),
		Rcpp::_["control"] = PGN_control,
		Rcpp::_["method"] = "BFGS",
		Rcpp::_["hessian"] = true,
		Rcpp::_["ind_List"] = ind_List,
		Rcpp::_["y_1"] = y_1,
		Rcpp::_["y_0"] = y_0,
		Rcpp::_["z_1"] = z_1,
		Rcpp::_["z_0"] = z_0,
		Rcpp::_["y_d_1"] = y_d_1,
		Rcpp::_["y_d_0"] = y_d_0,
		Rcpp::_["z_d_1"] = z_d_1,
		Rcpp::_["z_d_0"] = z_d_0,
		Rcpp::_["pol_degrees"] = pol_degrees,
		Rcpp::_["z_sd"] = z_sd,
		Rcpp::_["is_minus"] = true); //true because of minimization);

	//calculate additional values

		//get vector of estimated values

	NumericVector x1 = optim_results["par"];

	int x1_n = x1.size();

		//calcukate log-likelihood and AIC

	double lnL = optim_results["value"];
	lnL = lnL * (-1);

	double AIC = 2 * (x1_n - lnL);

		//get polynomial coefficients

	NumericVector pol_coefficients = NumericVector(pol_coefficients_n);

	if ((z_K != 0) | (y_K != 0))
	{
		pol_coefficients = x1[pol_coefficients_ind];
	}

	pol_coefficients.push_front(1);

		//get mean and sd values
	double z_mean = x1[z_mean_ind];
	double y_mean = x1[y_mean_ind];
	double y_sd = x1[y_sd_ind];

		//get coefficients

	NumericVector z_coef = x1[z_coef_ind];
	NumericVector y_coef = x1[y_coef_ind];

		//get covariance matrix

	NumericMatrix cov_mat = solve(optim_results["hessian"]);

	//Prepare beautifull results output

	NumericMatrix results(x1_n, 4);

	StringVector results_cols = StringVector::create("Estimate", "Std. Error", "z value", "P(>|z|)");

	StringVector results_rows(x1_n);

	//polIndex_mat

		//for alpha

	double z_stat = 0;
	NumericVector F_z_stat;
	for (int i = 1; i <= pol_coefficients_n; i++)
	{
		results_rows[(i - 1)] = "a_" + std::to_string((int)polIndex_mat(0, i)) + "_" + std::to_string((int)polIndex_mat(1, i));
		results((i - 1), 0) = pol_coefficients[(i - 1) + 1];
		results((i - 1), 1) = sqrt(cov_mat((i - 1), (i - 1)));
		z_stat = results((i - 1), 0) / results((i - 1), 1);
		F_z_stat = pnorm(NumericVector::create(z_stat));
		results((i - 1), 2) = z_stat;
		results((i - 1), 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}

		//for z_mean
		
	results_rows[z_mean_ind] = "z_mean";
	results(z_mean_ind, 0) = x1[z_mean_ind];
	results(z_mean_ind, 1) = sqrt(cov_mat((z_mean_ind - 1), (z_mean_ind - 1)));
	z_stat = results(z_mean_ind, 0) / results(z_mean_ind, 1);
	F_z_stat = pnorm(NumericVector::create(z_stat));
	results(z_mean_ind, 2) = z_stat;
	results(z_mean_ind, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);

		//for y_mean

	results_rows[y_mean_ind] = "y_mean";
	results(y_mean_ind, 0) = x1[y_mean_ind];
	results(y_mean_ind, 1) = sqrt(cov_mat((y_mean_ind), (y_mean_ind)));
	z_stat = results(y_mean_ind, 0) / results(y_mean_ind, 1);
	F_z_stat = pnorm(NumericVector::create(z_stat));
	results(y_mean_ind, 2) = z_stat;
	results(y_mean_ind, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);

		//for y_sd

	results_rows[y_sd_ind] = "y_sd";
	results(y_sd_ind, 0) = x1[y_sd_ind];
	results(y_sd_ind, 1) = sqrt(cov_mat((y_sd_ind), (y_sd_ind)));
	z_stat = results(y_sd_ind, 0) / results(y_sd_ind, 1);
	F_z_stat = pnorm(NumericVector::create(z_stat));
	results(y_sd_ind, 2) = z_stat;
	results(y_sd_ind, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);

		//for z coefficients

	for (int i = 0; i < (z_d_col - 1); i++)
	{
		results_rows[z_coef_ind[i]] = z_df_names(i + 2);
		results(z_coef_ind[i], 0) = x1[z_coef_ind[i]];
		results(z_coef_ind[i], 1) = sqrt(cov_mat(z_coef_ind[i], z_coef_ind[i]));
		z_stat = results(z_coef_ind[i], 0) / results(z_coef_ind[i], 1);
		F_z_stat = pnorm(NumericVector::create(z_stat));
		results(z_coef_ind[i], 2) = z_stat;
		results(z_coef_ind[i], 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}

		//for y coefficients

	for (int i = 0; i < y_d_col; i++)
	{
		results_rows[y_coef_ind[i]] = y_df_names(i + 1);
		results(y_coef_ind[i], 0) = x1[y_coef_ind[i]];
		results(y_coef_ind[i], 1) = sqrt(cov_mat(y_coef_ind[i], y_coef_ind[i]));
		double z_stat = results(y_coef_ind[i], 0) / results(y_coef_ind[i], 1);
		F_z_stat = pnorm(NumericVector::create(z_stat));
		results(y_coef_ind[i], 2) = z_stat;
		results(y_coef_ind[i], 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}

		//assign names to the output

	rownames(results) = results_rows;

	colnames(results) = results_cols;

	x1.names() = results_rows;

	rownames(cov_mat) = results_rows;
	colnames(cov_mat) = results_rows;

	if ((z_K != 0) & (y_K != 0))
	{
		pol_coefficients.names() = c_R("a_0", results_rows[pol_coefficients_ind]);
	}

	z_coef.names() = results_rows[z_coef_ind];

	y_coef.names() = results_rows[y_coef_ind];

	//Calculate expectation and variance
	NumericVector z_e = ehpa(NumericMatrix(1, 1), pol_coefficients, pol_degrees,
		LogicalVector{false, false}, LogicalVector{false, false},
		NumericVector::create(z_mean, y_mean), NumericVector::create(z_sd, y_sd),
		NumericVector::create(1, 0));

	NumericVector z_e_2 = ehpa(NumericMatrix(1, 1), pol_coefficients, pol_degrees,
		LogicalVector{false, false}, LogicalVector{false, false},
		NumericVector::create(z_mean, y_mean), NumericVector::create(z_sd, y_sd),
		NumericVector::create(2, 0));

	NumericVector y_e = ehpa(NumericMatrix(1, 1), pol_coefficients, pol_degrees,
		LogicalVector{false, false}, LogicalVector{false, false},
		NumericVector::create(z_mean, y_mean), NumericVector::create(z_sd, y_sd),
		NumericVector::create(0, 1));

	NumericVector y_e_2 = ehpa(NumericMatrix(1, 1), pol_coefficients, pol_degrees,
		LogicalVector{false, false}, LogicalVector{false, false},
		NumericVector::create(z_mean, y_mean), NumericVector::create(z_sd, y_sd),
		NumericVector::create(0, 2));

	NumericVector z_y_e = ehpa(NumericMatrix(1, 1), pol_coefficients, pol_degrees,
		LogicalVector{false, false}, LogicalVector{false, false},
		NumericVector::create(z_mean, y_mean), NumericVector::create(z_sd, y_sd),
		NumericVector::create(1, 1));

	double z_v = z_e_2[0] - z_e[0] * z_e[0];
	double y_v = y_e_2[0] - y_e[0] * y_e[0];
	double z_y_c = z_y_e[0] - z_e[0] * y_e[0];
	double rho = z_y_c / sqrt(z_v * y_v);

	List re_moments = List::create(Named("optim") = optim_results, 
		Named("selection_exp") = z_e,
		Named("outcome_exp") = y_e,
		Named("selection_var") = z_v,
		Named("outcome_var") = y_v,
		Named("errors_covariance") = z_y_c,
		Named("rho") = rho);

	List data_List = List::create(Named("data_z") = z_df,
		Named("data_y") = y_df,
		Named("dataframe") = data);

	List return_result = List::create(Named("optim") = optim_results,
		Named("x1") = x1,
		Named("Newey") = Newey,
		Named("log-likelihood") = lnL,
		Named("AIC") = AIC,
		Named("n_obs") = n,
		Named("data_List") = data_List,
		Named("results") = results,
		Named("z_mean") = z_mean,
		Named("y_mean") = y_mean, 
		Named("z_sd") = z_sd, 
		Named("y_sd") = y_sd,
		Named("pol_coefficients") = pol_coefficients,
		Named("pol_degrees") = pol_degrees,
		Named("y_coef") = x1[y_coef_ind],
		Named("z_coef") = x1[z_coef_ind],
		Named("selection_formula") = selection,
		Named("outcome_formula") = outcome,
		Named("re_moments") = re_moments,
		Named("ind_List") = ind_List);

	return_result.attr("class") = "hpaSelection";

	return(return_result);
}

// Perform log-likelihood function estimation for selection model
double hpaSelectionLnLOptim(NumericVector x0,
	List ind_List,
	arma::vec y_1,
	arma::vec y_0,
	arma::vec z_1,
	arma::vec z_0,
	arma::mat y_d_1,
	arma::mat y_d_0,
	arma::mat z_d_1,
	arma::mat z_d_0,
	NumericVector pol_degrees,
	double z_sd,
	bool is_minus = false)
{
	//Get values from the list

	NumericVector pol_coefficients_ind = ind_List["pol_coefficients_ind"];
	int z_mean_ind = ind_List["z_mean_ind"];
	int y_mean_ind = ind_List["y_mean_ind"];
	int y_sd_ind = ind_List["y_sd_ind"];
	NumericVector z_coef_ind = ind_List["z_coef_ind"];
	NumericVector y_coef_ind = ind_List["y_coef_ind"];
	

	//Assign estimated parameters values to corresponding vectors

		//polynomial coefficients and degrees
	NumericVector pol_coefficients = x0[pol_coefficients_ind];
	pol_coefficients.push_front(1); //add alpha(0...0)

		//mean
	double z_mean = x0[z_mean_ind];
	double y_mean = x0[y_mean_ind];
	NumericVector mean = {z_mean, y_mean};

		//sd
	double y_sd = x0[y_sd_ind];
	NumericVector sd = {z_sd, y_sd};

	if (y_sd <= 0)
	{
		return(999999999);
	}

		//coefficients for independend variables

	arma::vec y_coef = as<arma::vec>(x0[y_coef_ind]);
	NumericVector z_coef_temporal = x0[z_coef_ind];
	z_coef_temporal.push_front(1); //for fixed coefficient
	arma::vec z_coef = as<arma::vec>(z_coef_temporal);

		//get estimates for z*

	NumericVector z_h_1 = wrap(z_d_1 * z_coef);
	NumericVector z_h_0 = wrap(z_d_0 * z_coef);

		//get estimates for y and random errors

	arma::mat y_h_1 = y_d_1 * y_coef;

	NumericVector e_h_1 = wrap(y_1 - y_h_1);

	//concatenate e_h and z_h and prepare to insert into function

	NumericMatrix z_y_1 = NumericMatrix(z_h_1.size(), 2);
	NumericMatrix z_y_0 = NumericMatrix(z_h_0.size(), 2);

	z_y_1(_, 0) = (-1) * z_h_1;
	z_y_1(_, 1) = e_h_1;

	z_y_0(_, 0) = (-1) * z_h_0;

	//likelihood calculation

	double lnL_y_1 = 0;

	double lnL_z_y_1 = 0;

	double lnL_z_y_0 = 0;

	lnL_y_1 = (1 - 2 * is_minus) * sum(log(dhpa(z_y_1,
		pol_coefficients, pol_degrees,
		LogicalVector{false, false}, LogicalVector{true, false},
		mean, sd)));

	lnL_z_y_1 = (1 - 2 * is_minus) * sum(log(1 - phpa(z_y_1,
		pol_coefficients, pol_degrees,
		LogicalVector{false, true}, LogicalVector{false, false},
		mean, sd)));

	lnL_z_y_0 = (1 - 2 * is_minus) * sum(log(phpa(z_y_0,
		pol_coefficients, pol_degrees,
		LogicalVector{false, false}, LogicalVector{false, true},
		mean, sd)));

	return(lnL_y_1 + lnL_z_y_1 + lnL_z_y_0);
}

//' Predict outcome and selection equation values from hpaSelection model
//' @description This function predicts outcome and selection equation values from hpaSelection model.
//' @param object Object of class "hpaSelection"
//' @param method string value indicating prediction method based on hermite polynomial approximation "HPA" or Newey method "Newey".
//' @template newdata_Template
//' @param is_cond logical; if \code{TRUE} (default) then conditional predictions will be estimated. Otherwise unconditional predictions will be returned.
//' @param is_outcome logical; if \code{TRUE} (default) then predictions for selection equation will be estimated using "HPA" method.
//' Otherwise selection equation predictions (probabilities) will be returned.
//' @details Note that Newey method can't predict conditional outcomes for zero selection equation value. Conditional probabilities for selection equation
//' could be estimated only when dependent variable from outcome equation is observable.
//' @return This function returns the list which structure depends on \code{method}, \code{is_probit} and \code{is_outcome} values.
//' @export
// [[Rcpp::export]]
List predict_hpaSelection(List object, DataFrame newdata = R_NilValue, std::string method = "HPA", 
	bool is_cond = true, bool is_outcome = true)
{

	List model = object;

	//Add additional environments

	Rcpp::Environment stats_env("package:stats");
	Rcpp::Function model_frame = stats_env["model.frame"];
	Rcpp::Function na_omit_R = stats_env["na.omit"];
	Rcpp::Function na_pass = stats_env["na.pass"];

	Rcpp::Environment base_env("package:base");
	Rcpp::Function as_data_frame = base_env["as.data.frame"];

	//Working with Data

	List Newey = model["Newey"];

	double z_mean = model["z_mean"];
	double y_mean = model["y_mean"];
	double z_sd = model["z_sd"];
	double y_sd = model["y_sd"];

	NumericVector pol_degrees = model["pol_degrees"];
	NumericVector pol_coefficients = model["pol_coefficients"];

		//Check wheather new dataframe has been supplied

	Rcpp::Formula selection = model["selection_formula"];
	Rcpp::Formula outcome = model["outcome_formula"];

	DataFrame data;

	if (newdata.size() == 0)
	{
		List data_List = model["data_List"];
		newdata = as_data_frame(data_List["dataframe"]);
	}

	data = newdata;

	//Extract dataframe from formula

	DataFrame z_df = model_frame(Rcpp::_["formula"] = selection, Rcpp::_["data"] = newdata, Rcpp::_["na.action"] = na_pass);
	DataFrame y_df = model_frame(Rcpp::_["formula"] = outcome, Rcpp::_["data"] = newdata, Rcpp::_["na.action"] = na_pass);

	CharacterVector z_df_names = z_df.names();
	CharacterVector y_df_names = y_df.names();

	int z_df_n = z_df.size();
	int y_df_n = y_df.size();

	//Extract dependend variables

	NumericVector z = z_df[0]; //it is reference
	NumericVector y = y_df[0]; //it is reference

	int n = z.size();

	//Extract independend variable

	NumericMatrix z_d(n, z_df_n - 1);//-1 because there is no constant term
	NumericMatrix y_d(n, y_df_n - 1);//-1 because there is no constant term

	int z_d_col = z_d.ncol();
	int y_d_col = y_d.ncol();

	for (int i = 0; i < z_d_col; i++)
	{
		z_d(_, i) = NumericVector(z_df[i + 1]);
	}
	for (int i = 0; i < y_d_col; i++)
	{
		y_d(_, i) = NumericVector(y_df[i + 1]);
	}

	//calculate latent variables values

	NumericVector z_latent = NumericVector(n);

	NumericVector z_coef;
	
	if (method == "HPA")
	{
		z_coef = model["z_coef"];
	}

	if (method == "Newey")
	{
		z_coef = Newey["z_coef"];
	}

	for (int i = 0; i < z_d_col; i++)
	{
		z_latent = z_latent + z_d(_, i) * z_coef[i];
	}

	//calculate unconditional y values without constant
	
	NumericVector y_uncond = NumericVector(n);

	NumericVector y_coef;

	if (method == "HPA")
	{
		y_coef = model["y_coef"];
	}

	if (method == "Newey")
	{
		y_coef = Newey["y_coef"];
	}

	for (int i = 0; i < y_d_col; i++)
	{
		y_uncond = y_uncond + y_d(_, i) * y_coef[i];
	}

	//calculate conditional expectations

	List results;

	NumericVector y_cond = NumericVector(n);

	if (!is_outcome)
	{
		NumericMatrix e_mat = NumericMatrix(n, 1);
		e_mat(_, 0) = y - y_uncond;
		NumericMatrix z_y = NumericMatrix(n, 2);
		z_y(_, 0) = z_latent * (-1);
		z_y(_, 1) = e_mat;
		if (is_cond)
		{
			NumericVector z_prob = 1 - phpa(z_y, pol_coefficients,
				pol_degrees,
				LogicalVector::create(false, true), LogicalVector::create(false, false),
				NumericVector::create(z_mean, y_mean), NumericVector::create(z_sd, y_sd));
			return(List::create(Named("prob") = z_prob));
		}
		else
		{
			NumericVector z_prob = 1 - phpa(z_y, pol_coefficients,
				pol_degrees,
				LogicalVector::create(false, false), LogicalVector::create(false, true),
				NumericVector::create(z_mean, y_mean), NumericVector::create(z_sd, y_sd));
			return(List::create(Named("prob") = z_prob));
		}
	}

	if (method == "HPA")
	{
		if (!is_cond)
		{
			NumericVector errors_exp_vec = ehpa(NumericMatrix(1, 1),
				pol_coefficients, pol_degrees,
				LogicalVector::create(false, false),
				LogicalVector::create(false, false),
				NumericVector::create(z_mean, y_mean), 
				NumericVector::create(y_mean, y_sd),
				NumericVector::create(0, 1));
			double errors_exp = errors_exp_vec[0];

			return(List::create(Named("y") = y_uncond + errors_exp));
		}

		NumericVector NegInfVec = NumericVector(n);
		NumericVector PosInfVec = NumericVector(n);

		std::fill(NegInfVec.begin(), NegInfVec.end(), R_NegInf);
		std::fill(PosInfVec.begin(), PosInfVec.end(), R_PosInf);

		//for 1 outcome

		NumericMatrix lower_1 = NumericMatrix(n, 2);
		lower_1(_, 0) = (-1) * z_latent;
		lower_1(_, 1) = NegInfVec;

		NumericMatrix upper_1 = NumericMatrix(n, 2);
		upper_1(_, 0) = PosInfVec;
		upper_1(_, 1) = PosInfVec;

		NumericVector e_tr_1 = etrhpa(lower_1, upper_1,
			pol_coefficients, pol_degrees,
			NumericVector::create(z_mean, y_mean),
			NumericVector::create(z_sd, y_sd),
			NumericVector::create(0, 1));

		NumericVector y_cond_1 = y_uncond + e_tr_1;

		//for 0 outcome

		NumericMatrix lower_0 = NumericMatrix(n, 2);
		lower_0(_, 0) = NegInfVec;
		lower_0(_, 1) = NegInfVec;

		NumericMatrix upper_0 = NumericMatrix(n, 2);
		upper_0(_, 0) = (-1) * z_latent;
		upper_0(_, 1) = PosInfVec;

		NumericVector e_tr_0 = etrhpa(lower_0, upper_0,
			pol_coefficients, pol_degrees,
			NumericVector::create(z_mean, y_mean),
			NumericVector::create(z_sd, y_sd),
			NumericVector::create(0, 1));

		NumericVector y_cond_0 = y_uncond + e_tr_0;

		//aggregate result

		y_cond = y_cond_1;
		y_cond[z == 0] = y_cond_0[z == 0];

		results = List::create(Named("y") = y_cond,
			Named("y_1") = y_cond_1,
			Named("y_0") = y_cond_0);
	}

	if (method == "Newey")
	{

		if (!is_cond)
		{
			return(List::create(Named("y") = y_uncond));
		}

		y_cond = y_uncond;

		double z_exp = Newey["selection_exp"];
		double z_var = Newey["selection_var"];

		int pol_elements = Newey["pol_elements"];

		NumericVector y_coef_mills = Newey["inv_mills_coef"];
		z_latent = (z_latent - z_exp) / sqrt(z_var);
		NumericVector z_mills = dnorm(z_latent) / pnorm(z_latent);

		for (int i = 0; i <= pol_elements; i++)
		{
			y_cond = y_cond + pow(z_mills, i) * y_coef_mills[i];
		}

		results = List::create(Named("y_1") = y_cond);
	}

	return(results);
}

//' Summarizing hpaSelection Fits
//' @description This function summarizing hpaSelection Fits
//' @param object Object of class "hpaSelection"
//' @return This function returns the same list as \code{\link[hpa]{hpaSelection}} function changing it's class to "summary.hpaSelection".
//' @export
// [[Rcpp::export]]
List summary_hpaSelection(List object) 
{

	List return_result = clone(object); //in order to preserve model class

	return_result.attr("class") = "summary.hpaSelection";

	return(return_result);
}

//' Summary for hpaSelection output
//' @param x Object of class "hpaSelection"
//' @export
// [[Rcpp::export]]
void print_summary_hpaSelection(List x) {

	List model = x;

	//Load additional environments
	Rcpp::Environment base_env("package:base");
	Rcpp::Function as_table = base_env["as.table"];
	Rcpp::Function cbind = base_env["cbind"];
	Rcpp::Function round_R = base_env["round"];
	Rcpp::Function as_data_frame = base_env["as.data.frame"];

	NumericMatrix results = model["results"];
	results = round_R(Rcpp::_["x"] = results, Rcpp::_["digits"] = 5);

	//extract list of indexes
	List ind_List = model["ind_List"];
	NumericVector pol_coefficients_ind = ind_List["pol_coefficients_ind"];
	NumericVector z_mean_ind = ind_List["z_mean_ind"];
	NumericVector z_coef_ind = ind_List["z_coef_ind"];
	NumericVector y_mean_ind = ind_List["y_mean_ind"];
	NumericVector y_sd_ind = ind_List["y_sd_ind"];
	NumericVector y_coef_ind = ind_List["y_coef_ind"];

	//other stuff

	NumericVector x1 = model["x1"];

	List data_List = model["data_List"];

	DataFrame data = as_data_frame(data_List["dataframe"]);
	DataFrame data_z = as_data_frame(data_List["data_z"]);
	StringVector data_names_z = data_z.names();
	NumericVector z = data_z[0];
	int n_censored = sum(z);

	NumericVector p_values = results(_, 3);

	StringVector stars = starVector(p_values);

	NumericVector z_coef = model["z_coef"];
	NumericVector y_coef = model["y_coef"];

	StringVector results_rownames = rownames(results);
	StringVector results_colnames = colnames(results);

	double z_sd = model["z_sd"];

	double lnL = model["log-likelihood"];
	double AIC = model["AIC"];
	int n_obs = model["n_obs"];
	int df = x1.size();
	std::string lnL_string = "Log-Likelihood: " + std::to_string(lnL) + "\n";
	std::string AIC_string = "AIC: " + std::to_string(AIC) + "\n";
	std::string n_obs_string = "Observations: " + std::to_string(n_obs) + " (" + std::to_string(n_censored) + " observed)"+ "\n";
	std::string df_string = std::to_string(df) + " free parameters (df = " + std::to_string(n_obs - df) + ")" + "\n";

	Rprintf("%s", "--------------------------------------------------------------\n");

	Rprintf("%s", "Semi-nonparametric selection model estimation\n");

	std::printf("%s", lnL_string.c_str());
	std::printf("%s", AIC_string.c_str());
	std::printf("%s", n_obs_string.c_str());
	std::printf("%s", df_string.c_str());

	Rprintf("%s", "Selection equation coefficients:\n");
	int z_coef_first = z_coef_ind[0];
	int z_coef_last = z_coef_ind[z_coef_ind.size() - 1];
	NumericMatrix z_coef_results = results(Range(z_coef_first, z_coef_last), _);
	rownames(z_coef_results) = results_rownames[z_coef_ind];
	colnames(z_coef_results) = results_colnames;
	print(as_table(cbind(z_coef_results, stars[z_coef_ind])));

	Rprintf("%s", "---\n");

	Rprintf("%s", "Outcome equation coefficients:\n");
	int y_coef_first = y_coef_ind[0];
	int y_coef_last = y_coef_ind[y_coef_ind.size() - 1];
	NumericMatrix y_coef_results = results(Range(y_coef_first, y_coef_last), _);
	rownames(y_coef_results) = results_rownames[y_coef_ind];
	colnames(y_coef_results) = results_colnames;
	print(as_table(cbind(y_coef_results, stars[y_coef_ind])));

	Rprintf("%s", "---\n");

	Rprintf("%s", "Distribution parameters:\n");
	int distr_first = 0;
	int distr_last = df - z_coef_ind.size() - y_coef_ind.size() - 1;
	NumericMatrix distr_results = results(Range(distr_first, distr_last), _);
	StringVector distr_rownames = results_rownames[Range(distr_first, distr_last)];
	rownames(distr_results) = distr_rownames;
	colnames(distr_results) = results_colnames;
	print(as_table(cbind(distr_results, stars[Range(distr_first, distr_last)])));

	Rprintf("%s", "---\n");

	Rprintf("%s", "Selection equation fixed coefficients:\n");
	String new_str_names = data_names_z(1);
	std::string new_str_names_str = new_str_names;
	std::string new_str = new_str_names_str + " = 1" + "\n";
	Rprintf("%s", new_str.c_str());

	Rprintf("%s", "---\n");

	Rprintf("%s", "Fixed Distribution Parameters:\n");
	Rprintf("%s", "a_0 = 1\n");
	std::string new_str_z_sd = "z_sd = " + std::to_string(z_sd) + "\n";
	Rprintf("%s", new_str_z_sd.c_str());

	Rprintf("%s", "---\n");

	List re_moments = model["re_moments"];
	double rho = re_moments["rho"];
	std::string new_str_rho = "Correlation between random errors is " + std::to_string(rho) + "\n";
	Rprintf("%s", new_str_rho.c_str());


	Rprintf("%s", "---\n");

	Rprintf("%s", "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n");

	Rprintf("%s", "--------------------------------------------------------------\n");
}

//' Plot hpaSelection random errors approximated density
//' @param x Object of class "hpaSelection"
//' @param is_outcome logical; if TRUE then function plots the graph for outcome equation random errors. 
//' Otherwise plot for selection equation random errors will be plotted.
//' @return This function returns the list containing random error's expected value \code{errors_exp}
//' and variance \code{errors_var} estimates for selection (if \code{is_outcome = TRUE}) or outcome
//' (if \code{is_outcome = FALSE}) equation.
//' @export
// [[Rcpp::export]]
List plot_hpaSelection(List x, bool is_outcome = true) {

	List model = x;

	//Load additional environments
	Rcpp::Environment base_env("package:base");
	Rcpp::Function round_R = base_env["round"];

	Rcpp::Environment graphics_env("package:graphics");
	Rcpp::Function plot_R = graphics_env["plot"];

	//load data from the model

	double z_mean = model["z_mean"];
	double y_mean = model["y_mean"];
	double z_sd = model["z_sd"];
	double y_sd = model["y_sd"];

	NumericVector mean = NumericVector::create(z_mean, y_mean);
	NumericVector sd = NumericVector::create(z_sd, y_sd);

	NumericVector pol_degrees = model["pol_degrees"];
	NumericVector pol_coefficients = model["pol_coefficients"];

	//get random errors expectation and variance
	int eq_ind = is_outcome;

		NumericVector errors_exp_vec = ehpa(NumericMatrix(1, 1), 
			pol_coefficients, pol_degrees,
			LogicalVector::create(false, false), 
			LogicalVector::create(false, false),
			mean, sd,
			NumericVector::create(1 - eq_ind, eq_ind));
		double errors_exp = errors_exp_vec[0];

		NumericVector errors_exp_2_vec = ehpa(NumericMatrix(1, 1), 
			pol_coefficients, pol_degrees,
			LogicalVector::create(false, false), 
			LogicalVector::create(false, false),
			mean, sd,
			NumericVector::create(2 * (1 - eq_ind), 2 * eq_ind));
		double errors_exp_2 = errors_exp_2_vec[0];

		double errors_var = errors_exp_2 - (errors_exp * errors_exp);

	//adjust precision

	double plot_min = errors_exp - 3 * sqrt(errors_var);
	double plot_max = errors_exp + 3 * sqrt(errors_var);

	int n = 10000;

	double precise = (plot_max - plot_min) / n;

	NumericMatrix x_matr = NumericMatrix(n, 2);
	x_matr(0, 0) = plot_min;
	x_matr(0, 1) = plot_min;

	for (int i = 1; i < n; i++)
	{
		x_matr(i, eq_ind) = x_matr(i - 1, eq_ind) + precise;
	}

	//calculate densities

	NumericVector den = dhpa(x_matr,
		pol_coefficients, pol_degrees,
		LogicalVector::create(false, false),
		LogicalVector::create(is_outcome, !is_outcome),
		mean, sd);

	double den_min = min(den);
	double den_max = max(den);

	NumericVector x_vec = x_matr(_, eq_ind);

	plot_R(Rcpp::_["x"] = x_vec, Rcpp::_["y"] = den,
		Rcpp::_["xlim"] = NumericVector::create(plot_min, plot_max),
		Rcpp::_["xaxp"] = NumericVector::create(plot_min, plot_max, 5),
		Rcpp::_["yaxp"] = NumericVector::create(den_min, den_max, 5),
		Rcpp::_["type"] = "l",
		Rcpp::_["lwd"] = 3,
		Rcpp::_["main"] = "Random Errors Density Approximation Plot",
		Rcpp::_["xlab"] = "value",
		Rcpp::_["ylab"] = "density");

	List moments = List::create(Named("errors_exp") = errors_exp,
		Named("errors_var") = errors_var);

	return(moments);
}

//' Calculates AIC for "hpaSelection" object
//' @description This function calculates AIC for "hpaSelection" object
//' @param object Object of class "hpaSelection"
//' @template AIC_Template
//' @export
// [[Rcpp::export]]
double AIC_hpaSelection(List object, double k = 2) 
{

	double AIC = object["AIC"];

	if (k == 2)
	{
		return(AIC);
	}

	NumericVector x1 = object["x1"];

	AIC += (k - 2) * x1.size();

	return(AIC);
}

//' Calculates log-likelihood for "hpaSelection" object
//' @description This function calculates log-likelihood for "hpaSelection" object
//' @param object Object of class "hpaSelection"
//' @export
// [[Rcpp::export]]
double logLik_hpaSelection(List object) 
{

	double lnL = object["log-likelihood"];

	return(lnL);
}
