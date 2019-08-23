#include "hpaMain.h"
#include "hpaML.h"
#include "hpaBinary.h"
#include "polynomialIndex.h"

#include <RcppArmadillo.h>
using namespace RcppArmadillo;

// [[Rcpp::depends(RcppArmadillo)]]

//' Perform semi-nonparametric binary choice model estimation
//' @description This function performs semi-nonparametric single index binary choice model estimation
//' via hermite polynomial densities approximation.
//' @template formula_Template
//' @template data_Template
//' @template K_Template
//' @template z_mean_fixed_Template
//' @template z_sd_fixed_Template
//' @template z_constant_fixed_Template
//' @template z_coef_first_fixed_Template
//' @template x0_binary_Template
//' @param is_x0_probit logical; if \code{TRUE} (default) then initial points for optimization routine will be
//' obtained by probit model estimated via \link[stats]{glm} function.
//' @template is_sequence_Template
//' @template hpa_likelihood_details_Template
//' @template GN_details_Template
//' @template first_coef_Template
//' @details Note that if \code{is_z_coef_first_fixed} value is TRUE then the coefficient for the first
//' independent variable in \code{formula} will be fixed to 1.
//' @template sd_adjust_Template
//' @template is_numeric_Template
//' @template parametric_paradigm_Template
//' @template optim_details_Template
//' @return This function returns an object of class "hpaBinary".\cr \cr
//' An object of class "hpaBinary" is a list containing the following components:
//' \itemize{
//' \item \code{optim} - \code{\link[stats]{optim}} function output.
//' \item \code{x1} - numeric vector of distribution parameters estimates.
//' \item \code{mean} - mean (mu) parameter of density function estimate.
//' \item \code{sd} - sd (sigma) parameter of density function estimate.
//' \item \code{pol_coefficients} - polynomial coefficients estimates.
//' \item \code{pol_degrees} - the same as \code{K} input parameter.
//' \item \code{coefficients} - regression (single index) coefficients estimates.
//' \item \code{cov_matrix} - estimated parameters covariance matrix estimate.
//' \item \code{results} - numeric matrix representing estimation results.
//' \item \code{log-likelihood} - value of Log-Likelihood function.
//' \item \code{AIC} - AIC value.
//' \item \code{errors_exp} - random error expectation estimate.
//' \item \code{errors_var} - random error variance estimate.
//' \item \code{dataframe} - dataframe containing variables mentioned in \code{formula} without \code{NA} values.
//' \item \code{model_Lists} - lists containing information about fixed parameters and parameters indexes in \code{x1}.
//' \item \code{n_obs} - number of observations.
//' \item \code{z_latent} - latent variable (signle index) estimates.
//' \item \code{z_prob} - probabilities of positive outcome (i.e. 1) estimates.}
//' @seealso \link[hpa]{summary.hpaBinary}, \link[hpa]{predict.hpaBinary}, \link[hpa]{plot.hpaBinary},
//' \link[hpa]{AIC.hpaBinary}, \link[hpa]{logLik.hpaBinary}
//' @template hpaBinary_examples_Template
//' @export	
// [[Rcpp::export]]
List hpaBinary(Rcpp::Formula formula,
	DataFrame data,
	int K = 1,
	double z_mean_fixed = NA_REAL,
	double z_sd_fixed = NA_REAL,
	double z_constant_fixed = 0,
	bool is_z_coef_first_fixed = true,
	bool is_x0_probit = true,
	bool is_sequence = false,
	NumericVector x0 = NumericVector(0)) 
{

	//Load additional environments

	Rcpp::Environment stats_env("package:stats");
	Rcpp::Function optim = stats_env["optim"];
	Rcpp::Function model_frame = stats_env["model.frame"];
	Rcpp::Function na_omit_R = stats_env["na.omit"];
	Rcpp::Function glm = stats_env["glm"];
	Rcpp::Function binomial = stats_env["binomial"];
	Rcpp::Function na_pass = stats_env["na.pass"];

	Rcpp::Environment base_env("package:base");
	Rcpp::Function solve = base_env["solve"];
	Rcpp::Function class_R = base_env["class"];
	Rcpp::Function c_R = base_env["c"];

	//Initialize polynomial structure related values

	int pol_coefficients_n = K;

	//Initialize bool values related to fixed parameters

	bool is_z_mean_fixed = !R_IsNA(z_mean_fixed);

	bool is_z_sd_fixed = !R_IsNA(z_sd_fixed);

	bool is_z_constant_fixed = !R_IsNA(z_constant_fixed);

	//Working with Data

		//Extract dataframe from formula

	DataFrame z_df = model_frame(Rcpp::_["formula"] = formula, Rcpp::_["data"] = data, Rcpp::_["na.action"] = na_pass);

	z_df = na_omit_R(z_df);

	int z_df_n = z_df.size();

		//Extract binary dependend variable
	 
	NumericVector z = z_df[0]; //it is reference

	int n = z.size();

		//Extract independend variable

	NumericMatrix z_d(n, z_df_n - is_z_constant_fixed);

	int z_d_col = z_d.ncol();

	if (!is_z_constant_fixed)
	{
		z_d(_, 0) = (NumericVector(n) + 1); //add constant

	}

	for (int i = 0; i < (z_d_col - !is_z_constant_fixed); i++)
	{
		z_d(_, i + !is_z_constant_fixed) = as<NumericVector>(z_df[i + 1]);
	}

	//Sequential estimation
	if (is_sequence)
	{
		//for K=0
		List results(K + 1);
		results[0] = hpaBinary(formula, data, 0, z_mean_fixed, z_sd_fixed, z_constant_fixed, is_z_coef_first_fixed, true, false);
		List results_current = results[0];
		NumericVector x1 = results_current["x1"];
		int x0_n = x1.size() + 1; //add one more alpha parameter for the next estimation
		NumericVector x0 = NumericVector(x0_n);
		for (int i = 1; i < x0_n; i++)
		{
			x0[i] = x1[i - 1];
		}
		//for other K
		for (int i = 1; i <= K; i++)
		{
			if (is_z_sd_fixed)
			{
				z_sd_fixed = results_current["sd"];
			}
			results[i] = hpaBinary(formula, data, i, z_mean_fixed, z_sd_fixed, z_constant_fixed, is_z_coef_first_fixed, false, false, x0);
			results_current = results[i];
			x1 = results_current["x1"];
			x0_n++;
			x0 = NumericVector(x0_n);
			for (int j = 0; j < x0_n; j++)
			{
				if (i > j)
				{
					x0[j] = x1[j];
				}
				if (i < j)
				{
					x0[j] = x1[j - 1];
				}
			}
		}
		return(results);
	}

	//Create initial values vector
	bool x0_given = true;

	if (x0.size() == 0)
	{
		x0_given = false;
		//x0 dimensions are equal to the estimated (nonfixed) parameters number
		//note thet constant is already in z_d_col or not
		x0 = NumericVector(pol_coefficients_n +
							!is_z_mean_fixed + !is_z_sd_fixed +
							z_d_col - is_z_coef_first_fixed);
	}

	//Assign indexes

		//Initialize additional index and upper value for some loops

	int k = 0; //to account for fixed parameters

	int lower_ind = 0;

	int upper_ind = pol_coefficients_n;

		//for polynomial coefficients

	NumericVector pol_coefficients_ind(pol_coefficients_n);

	if (K != 0)
	{
		for (int i = lower_ind; i < upper_ind; i++)
		{
			pol_coefficients_ind[i] = i;
		}
	} else {
		pol_coefficients_ind = NumericVector::create(0);
	}
		//for mean vector

	int z_mean_ind = 0;

	if (!is_z_mean_fixed)
	{
		z_mean_ind = pol_coefficients_n + k;
		k++;
	}

		//for sd vector
	int z_sd_ind = 0;

	if (!is_z_sd_fixed)
	{
		z_sd_ind = pol_coefficients_n + k;
		k++;
	}
		//for z coefficients
		//note that z_d_col may contain or not the constant term

	NumericVector z_coef_ind(z_d_col - is_z_coef_first_fixed);

	lower_ind = pol_coefficients_n + k;

	upper_ind = pol_coefficients_n + z_coef_ind.size() + k;

	for (int i = lower_ind; i < upper_ind; i++)
	{
		z_coef_ind[i - lower_ind] = i;
	}

	//Convert to arma

	arma::vec z_arma = as<arma::vec>(z);

	arma::mat z_d_arma = as<arma::mat>(z_d);

	//Divide into 0 and 1 samples

	arma::vec z_1 = as<arma::vec>(z[z == 1]);

	arma::mat z_d_1 = (as<arma::mat>(z_d)).rows(arma::find(z_arma == 1));

	arma::vec z_0 = as<arma::vec>(z[z == 0]);

	arma::mat z_d_0 = (as<arma::mat>(z_d)).rows(arma::find(z_arma == 0));

	//set initial sd value to 1

	if (!is_z_sd_fixed & !x0_given)
	{
		x0[z_sd_ind] = 1;
	}

	//Estimate intial values via probit model using glm function

	int z_coef_n = z_coef_ind.size();

	if (is_x0_probit & !x0_given)
	{
		List model_probit = glm(Rcpp::_["formula"] = formula, Rcpp::_["data"] = data,
			Rcpp::_["family"] = binomial(Rcpp::_["link"] = "probit"));
		NumericVector glm_coef = model_probit["coefficients"];
		double coef_adjust = std::abs(glm_coef[1]);
		if (is_z_coef_first_fixed)
		{
			//addjust sd to first coefficient value
			if (is_z_sd_fixed)
			{
				z_sd_fixed = z_sd_fixed / coef_adjust;
			}
			else
			{
				x0[z_sd_ind] = x0[z_sd_ind] / coef_adjust;
			}
			glm_coef = glm_coef / coef_adjust;
		}
		else
		{
			if (is_z_sd_fixed)
			{
				glm_coef = glm_coef / z_sd_fixed;
			}
		}
		if (!is_z_constant_fixed)
		{
			x0[z_coef_ind[0]] = glm_coef[0]; //already adjusted because glm_coef = glm_coef / coef_adjust 
		}
		else
		{
			if (!is_z_mean_fixed)
			{
				x0[z_mean_ind] = glm_coef[0]; //already adjusted because glm_coef = glm_coef / coef_adjust 
			}
			else
			{
				z_mean_fixed = glm_coef[0];
			}
		}
		for (int i = 0; i < z_coef_n; i++)
		{
			x0[z_coef_ind[i+!is_z_constant_fixed]] = glm_coef[i + is_z_constant_fixed + is_z_coef_first_fixed];
		}
	}

	//Create list for some variables because unfortunatelly optim function has limitation for the
	//parameters number (parameters itself not estimated)

	List is_List = List::create(Named("is_z_coef_first_fixed") = is_z_coef_first_fixed, 
		Named("is_z_mean_fixed") = is_z_mean_fixed,
		Named("is_z_sd_fixed") = is_z_sd_fixed,
		Named("is_z_constant_fixed") = is_z_constant_fixed
	);

	//Apply optimization routine

	List PGN_control = List::create(Named("maxit") = 1000000);

	List optim_results = optim(
		Rcpp::_["par"] = x0,
		Rcpp::_["fn"] = Rcpp::InternalFunction(&hpaBinaryLnLOptim),
		Rcpp::_["control"] = PGN_control,
		Rcpp::_["method"] = "BFGS",
		Rcpp::_["hessian"] = true,
		Rcpp::_["is_List"] = is_List,
		Rcpp::_["z_1"] = z_1,
		Rcpp::_["z_0"] = z_0,
		Rcpp::_["z_d_1"] = z_d_1,
		Rcpp::_["z_d_0"] = z_d_0,
		Rcpp::_["K"] = K,
		Rcpp::_["z_mean_fixed"] = z_mean_fixed,
		Rcpp::_["z_sd_fixed"] = z_sd_fixed,
		Rcpp::_["z_constant_fixed"] = z_constant_fixed,
		Rcpp::_["pol_coefficients_ind"] = pol_coefficients_ind,
		Rcpp::_["z_mean_ind"] = z_mean_ind,
		Rcpp::_["z_sd_ind"] = z_sd_ind,
		Rcpp::_["z_coef_ind"] = z_coef_ind,
		Rcpp::_["is_minus"] = true); //true because of minimization);

	//extract coefficients and covariance matrix

	NumericVector x1 = optim_results["par"];

	double optim_value = optim_results["value"];

	int x1_n = x0.size();

	NumericVector pol_coefficients = NumericVector(K);

	if (K != 0)
	{
		pol_coefficients = x1[pol_coefficients_ind];
	}

	pol_coefficients.push_front(1);

	NumericVector z_mean = NumericVector(1);
	if(!is_z_mean_fixed)
	{ 
		z_mean = NumericVector::create(x1[z_mean_ind]);
	}
	else
	{
		z_mean = NumericVector::create(z_mean_fixed);
	}

	NumericVector z_sd = NumericVector(1);
	if (!is_z_sd_fixed)
	{
		z_sd = NumericVector::create(x1[z_sd_ind]);
	}
	else
	{
		z_sd = NumericVector::create(z_sd_fixed);
	}

	NumericVector z_coef = x1[z_coef_ind];

	NumericMatrix cov_mat = solve(optim_results["hessian"]);

	//Prepare beautifull results output

	NumericMatrix results(x1_n, 4);

	StringVector results_cols = StringVector::create("Estimate", "Std. Error", "z value", "P(>|z|)");

	StringVector results_rows(x1_n);

	NumericMatrix pol_ind = polynomialIndex(NumericVector::create(K));

		//for alpha

		for (int i = 1; i < (K+1); i++)
		{
			//Convert double to string
			std::stringstream ss;
			ss << i;
			std::string str2 = ss.str();
			//
			results_rows[i - 1] = "a_" + str2;
			results((i - 1), 0) = pol_coefficients[i];
			results((i - 1), 1) = sqrt(cov_mat((i - 1), (i - 1)));
			double z_stat = results((i - 1), 0) / results((i - 1), 1);
			NumericVector F_z_stat = pnorm(NumericVector::create(z_stat));
			results((i - 1), 2) = z_stat;
			results((i - 1), 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
		}

		//for mean
	if (!is_z_mean_fixed)
	{
		results_rows[z_mean_ind] = "mean";
		results(z_mean_ind, 0) = x1[z_mean_ind];
		results(z_mean_ind, 1) = sqrt(cov_mat((z_mean_ind), (z_mean_ind)));
		double z_stat = results(z_mean_ind, 0) / results(z_mean_ind, 1);
		NumericVector F_z_stat = pnorm(NumericVector::create(z_stat));
		results(z_mean_ind, 2) = z_stat;
		results(z_mean_ind, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}

		//for sd
	if (!is_z_sd_fixed)
	{
		results_rows[z_sd_ind] = "sd";
		results(z_sd_ind, 0) = x1[z_sd_ind];
		results(z_sd_ind, 1) = sqrt(cov_mat((z_sd_ind), (z_sd_ind)));
		double z_stat = results(z_sd_ind, 0) / results(z_sd_ind, 1);
		NumericVector F_z_stat = pnorm(NumericVector::create(z_stat));
		results(z_sd_ind, 2) = z_stat;
		results(z_sd_ind, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}
	
		//for z coefficients
	CharacterVector z_df_names = z_df.names();

	if (!is_z_constant_fixed)
	{
		int z_constant_ind = K + !is_z_mean_fixed + !is_z_sd_fixed;
		results_rows(z_constant_ind) = "(Intercept)";
		results(z_constant_ind, 0) = z_coef[0];
		results(z_constant_ind, 1) = sqrt(cov_mat(z_coef_ind[0], z_coef_ind[0]));
		double z_stat = results(z_constant_ind, 0) / results(z_constant_ind, 1);
		NumericVector F_z_stat = pnorm(NumericVector::create(z_stat));
		results(z_constant_ind, 2) = z_stat;
		results(z_constant_ind, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}

	if (!is_z_coef_first_fixed)
	{
		int z_first_ind = K + !is_z_mean_fixed + !is_z_sd_fixed + !is_z_constant_fixed;
		results_rows(z_first_ind) = z_df_names[1];
		results(z_first_ind, 0) = z_coef[!is_z_constant_fixed];
		results(z_first_ind, 1) = sqrt(cov_mat(z_coef_ind[!is_z_constant_fixed], z_coef_ind[!is_z_constant_fixed]));
		double z_stat = results(z_first_ind, 0) / results(z_first_ind, 1);
		NumericVector F_z_stat = pnorm(NumericVector::create(z_stat));
		results(z_first_ind, 2) = z_stat;
		results(z_first_ind, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}


	int z_rows_k = K + !is_z_mean_fixed + !is_z_sd_fixed + !is_z_constant_fixed + !is_z_coef_first_fixed;
	int z_names_k = !is_z_constant_fixed + !is_z_coef_first_fixed; 
	for (int i = z_rows_k; i < x1_n; i++)
	{
		results_rows(i) = z_df_names[(z_names_k - !is_z_constant_fixed + is_z_coef_first_fixed) + 1];//+1 because the first is dependend variable
		results(i, 0) = z_coef[z_names_k];
		results(i, 1) = sqrt(cov_mat(z_coef_ind[z_names_k], z_coef_ind[z_names_k]));
		double z_stat = results(i, 0) / results(i, 1);
		NumericVector F_z_stat = pnorm(NumericVector::create(z_stat));
		results(i, 2) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
		results(i, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
		z_names_k++;
	}

		//for expectation and variance
	NumericVector z_e = ehpa(NumericMatrix(1, 1), pol_coefficients, NumericVector::create(K),
		LogicalVector::create(false), LogicalVector::create(false),
		z_mean, z_sd, 
		NumericVector::create(1));

	NumericVector z_e_2 = ehpa(NumericMatrix(1, 1), pol_coefficients, NumericVector::create(K),
		LogicalVector::create(0), LogicalVector::create(0),
		z_mean, z_sd,
		NumericVector::create(2));

		//assign names to the output

	rownames(results) = results_rows;

	colnames(results) = results_cols;

	x1.names() = results_rows;

	rownames(cov_mat) = results_rows;
	colnames(cov_mat) = results_rows;

	if (K != 0)
	{
		pol_coefficients.names() = c_R("a_0", results_rows[pol_coefficients_ind]);
	}

	z_coef.names() = results_rows[z_coef_ind];

	//Estimate latent variable and probabilities

		//coefficients for independend variables

	if (is_z_coef_first_fixed)
	{
		z_coef.push_front(1); //add identity coefficient for fixed value
		if (!is_z_constant_fixed) //change places of constant and first z_d coefficient
		{
			double z_coef_0 = z_coef[0];
			z_coef[0] = z_coef[1];
			z_coef[1] = z_coef_0;
		}
	}

	arma::vec z_coef_arma = as<arma::vec>(z_coef);

		//get estimates for z*

	NumericMatrix z_latent = wrap(z_d_arma * z_coef_arma);

	if (is_z_constant_fixed)
	{
		z_latent = z_latent + z_constant_fixed;
	}

	NumericVector z_prob = 1 - phpa((-1) * z_latent, pol_coefficients,
		NumericVector::create(K),
		LogicalVector::create(0), LogicalVector::create(0),
		z_mean, z_sd);

	//return results
	
	List ind_List = List::create(Named("pol_coefficients_ind") = pol_coefficients_ind,
		Named("z_mean_ind") = z_mean_ind,
		Named("z_sd_ind") = z_mean_ind,
		Named("z_coef_ind") = z_coef_ind);

	List fixed_List = List::create(
		Named("z_mean_fixed") = z_mean_fixed,
		Named("z_sd_fixed") = z_sd_fixed,
		Named("z_constant_fixed") = z_constant_fixed);

	List model_Lists = List::create(Named("is_List") = is_List,
		Named("ind_List") = ind_List,
		Named("fixed_List") = fixed_List);

	List return_result = List::create(Named("optim") = optim_results,
		Named("x1") = x1,
		Named("mean") = z_mean,
		Named("sd") = z_sd,
		Named("pol_coefficients") = pol_coefficients,
		Named("pol_degrees") = NumericVector::create(K),
		Named("coefficients") = z_coef,
		Named("results") = results,
		Named("errors_exp") = z_e[0],
		Named("errors_var") = (z_e_2[0] - z_e[0] * z_e[0]),
		Named("log-likelihood") = (-1) * optim_value,
		Named("AIC") = 2 * (x1_n + optim_value),
		Named("n_obs") = z_latent.nrow(),
		Named("z_latent") = z_latent,
		Named("z_prob") = z_prob,
		Named("formula") = formula,
		Named("dataframe") = z_df,
		Named("model_Lists") = model_Lists,
		Named("cov_mat") = cov_mat);

	return_result.attr("class") = "hpaBinary";

		//estimate probabilities 

	return(return_result);
}

// Perform semi-nonparametric log-likelihood function estimation for binary choice model
double hpaBinaryLnLOptim(NumericVector x0,
	List is_List,
	arma::vec z_1,
	arma::vec z_0,
	arma::mat z_d_1,
	arma::mat z_d_0,
	int K = 1,
	double z_mean_fixed = NA_REAL,
	double z_sd_fixed = NA_REAL,
	double z_constant_fixed = 0,
	NumericVector pol_coefficients_ind = NumericVector(0),
	int z_mean_ind = 1,
	int z_sd_ind = 2,
	NumericVector z_coef_ind = NumericVector(0),
	bool is_minus = false) {

	//Get values from the is_List

	bool is_z_coef_first_fixed = is_List["is_z_coef_first_fixed"];

	bool is_z_mean_fixed = is_List["is_z_mean_fixed"];

	bool is_z_sd_fixed = is_List["is_z_sd_fixed"];

	bool is_z_constant_fixed = is_List["is_z_constant_fixed"];

	//Assign estimated parameters values to corresponding vectors

		//polynomial coefficients and degrees

	NumericVector pol_coefficients = NumericVector(K);

	if (K != 0)
	{
		pol_coefficients = x0[pol_coefficients_ind];
	} 

	pol_coefficients.push_front(1);

	NumericVector pol_degrees = NumericVector(1);

	pol_degrees[0] = K;

		//mean value

	NumericVector z_mean = NumericVector(1);

	if (is_z_mean_fixed)
	{
		z_mean[0] = z_mean_fixed;
	}
	else {
		z_mean[0] = x0[z_mean_ind];
	}

		//sd value

	NumericVector z_sd = NumericVector(1);

	if (is_z_sd_fixed)
	{
		z_sd[0] = z_sd_fixed;
	}
	else {
		z_sd[0] = x0[z_sd_ind];
	}

		//coefficients for independend variables

	NumericVector z_coef_R = x0[z_coef_ind];

	if (is_z_coef_first_fixed)
	{
		z_coef_R.push_front(1); //add identity coefficient for fixed value
		if (!is_z_constant_fixed) //change places of constant and first z_d coefficient
		{
			double z_coef_0 = z_coef_R[0];
			z_coef_R[0] = z_coef_R[1];
			z_coef_R[1] = z_coef_0;
		}
	}

	arma::vec z_coef = as<arma::vec>(z_coef_R);

	//get estimates for z*

	NumericMatrix z_h_1 = wrap(-z_d_1 * z_coef);

	if (is_z_constant_fixed)
	{
		z_h_1 = z_h_1 - z_constant_fixed;
	}

	NumericMatrix z_h_0 = wrap(-z_d_0 * z_coef);

	if (is_z_constant_fixed)
	{
		z_h_0 = z_h_0 - z_constant_fixed;
	}

	//likelihood calculation

	double lnL_z_1 = 0;

	double lnL_z_0 = 0;

	lnL_z_1 = (1 - 2 * is_minus) * (sum(log(1 - phpa(z_h_1,
		pol_coefficients, pol_degrees,
		LogicalVector{false}, LogicalVector{false},
		z_mean, z_sd))));

	lnL_z_0 = (1 - 2 * is_minus) * sum(log(phpa(z_h_0,
		pol_coefficients, pol_degrees,
		LogicalVector{false}, LogicalVector{false},
		z_mean, z_sd)));

	double result_value = lnL_z_1 + lnL_z_0;

	return(result_value);
}

//' Predict method for hpaBinary
//' @param object Object of class "hpaBinary"
//' @template newdata_Template
//' @param is_prob logical; if TRUE (default) then function returns predicted probabilities. Otherwise latent variable
//' (single index) estimates will be returned.
//' @return This function returns predicted probabilities based on \code{\link[hpa]{hpaBinary}} estimation results.
//' @export
// [[Rcpp::export]]
NumericVector predict_hpaBinary(List object, DataFrame newdata = R_NilValue, bool is_prob = true)
{

	List model = object;

	//Add additional environments

	Rcpp::Environment stats_env("package:stats");
	Rcpp::Function model_frame = stats_env["model.frame"];
	Rcpp::Function na_omit_R = stats_env["na.omit"];

	Rcpp::Environment base_env("package:base");
	Rcpp::Function as_data_frame = base_env["as.data.frame"];

	//Extract variables from model

		//extract is values

	List model_Lists = model["model_Lists"];

	List is_List = model_Lists["is_List"];

	bool is_z_constant_fixed = is_List["is_z_constant_fixed"];

		//extract fixed values

	List fixed_List = model_Lists["fixed_List"];

	double z_constant_fixed = fixed_List["z_constant_fixed"];

		//extract coefficients

	NumericVector pol_coefficients = model["pol_coefficients"];

	NumericVector z_mean = model["mean"];

	NumericVector z_sd = model["sd"];

	NumericVector z_coef = model["coefficients"];

		//extract polynomial coefficients

	double K = pol_coefficients.size() - 1;

	//Check wheather new dataframe has been supplied

	DataFrame data = newdata;

	if (newdata.size() == 0)
	{
		newdata = as_data_frame(model["dataframe"]);
	}

	//Remove NA values

	data = na_omit_R(newdata);

	//Working with Data

		//Extract dataframe from formula

	Formula formula = model["formula"];

	DataFrame z_df = model_frame(Rcpp::_["formula"] = formula, Rcpp::_["data"] = data);

	int z_df_n = z_df.size();

	//Extract binary dependend variable

	NumericVector z = z_df[0]; //it is reference

	int n = z.size();

	//Extract independend variable

	NumericMatrix z_d(n, z_df_n - is_z_constant_fixed);

	if (!is_z_constant_fixed)
	{
		z_d(_, 0) = (NumericVector(n) + 1); //add constant
	}

	for (int i = 1; i < z_df_n; i++)
	{
		z_d(_, i - is_z_constant_fixed) = as<NumericVector>(z_df[i]);
	}

		//Convert to arma

	arma::vec z_arma = as<arma::vec>(z);

	arma::mat z_d_arma = as<arma::mat>(z_d);

	//Estimate latent variable and probabilities

		//coefficients for independend variables

	arma::vec z_coef_arma = as<arma::vec>(z_coef);

	//get estimates for z*

	NumericMatrix z_latent = wrap(z_d_arma * z_coef_arma);

	if (is_z_constant_fixed)
	{
		z_latent = z_latent + z_constant_fixed;
	}

	if (!is_prob)
	{
		NumericVector z_latent_vec = z_latent(_, 0);
		return(z_latent_vec);
	}

	NumericVector z_prob = 1 - phpa((-1) * z_latent, pol_coefficients,
		NumericVector::create(K),
		LogicalVector::create(0), LogicalVector::create(0),
		z_mean, z_sd);

	return(z_prob);
}

//' Summarizing hpaBinary Fits
//' @param object Object of class "hpaBinary"
//' @return This function returns the same list as \code{\link[hpa]{hpaBinary}} function changing it's class to "summary.hpaBinary".
//' @export
// [[Rcpp::export]]
List summary_hpaBinary(List object)
{

	List return_result = clone(object); //in order to preserve model class

	return_result.attr("class") = "summary.hpaBinary";

	return(return_result);
}

//' Summary for hpaBinary output
//' @param x Object of class "hpaML"
//' @export	
// [[Rcpp::export]]
void print_summary_hpaBinary(List x)
{

	//Extract the model

	List model = x;

	//Load additional environments
	Rcpp::Environment base_env("package:base");
	Rcpp::Function as_table = base_env["as.table"];
	Rcpp::Function cbind = base_env["cbind"];
	Rcpp::Function round_R = base_env["round"];
	Rcpp::Function as_data_frame = base_env["as.data.frame"];

	NumericMatrix results = model["results"];
	results = round_R(Rcpp::_["x"] = results, Rcpp::_["digits"] = 5);

	List model_Lists = model["model_Lists"];

	List ind_List = model_Lists["ind_List"];
	NumericVector z_coef_ind = ind_List["z_coef_ind"];

	//Extract is values

	List is_List = model_Lists["is_List"];
	bool is_z_coef_first_fixed = is_List["is_z_coef_first_fixed"];
	bool is_z_mean_fixed = is_List["is_z_mean_fixed"];
	bool is_z_sd_fixed = is_List["is_z_sd_fixed"];
	bool is_z_constant_fixed = is_List["is_z_constant_fixed"];

	//Extract fixed values

	List fixed_List = model_Lists["fixed_List"];
	double z_constant_fixed = fixed_List["z_constant_fixed"];

	//Other stuff

	NumericVector x1 = model["x1"];

	DataFrame data = as_data_frame(model["dataframe"]);
	StringVector data_names = data.names();

	NumericVector p_values = results(_, 3);

	StringVector stars = starVector(p_values);

	NumericVector z_coef = model["coefficients"];

	StringVector results_rownames = rownames(results);
	StringVector results_colnames = colnames(results);

	double mean = model["mean"];
	double sd = model["sd"];

	double lnL = model["log-likelihood"];
	double AIC = model["AIC"];
	int n_obs = model["n_obs"];
	int df = x1.size();
	std::string lnL_string = "Log-Likelihood: " + std::to_string(lnL) + "\n";
	std::string AIC_string = "AIC: " + std::to_string(AIC) + "\n";
	std::string n_obs_string = "Observations: " + std::to_string(n_obs) + "\n";
	std::string df_string = std::to_string(df) + " free parameters (df = " + std::to_string(n_obs - df) + ")" + "\n";

	Rprintf("%s", "--------------------------------------------------------------\n");

	Rprintf("%s", "Semi-nonparametric binary choice model estimation\n");

	std::printf("%s", lnL_string.c_str());
	std::printf("%s", AIC_string.c_str());
	std::printf("%s", n_obs_string.c_str());
	std::printf("%s", df_string.c_str());

	Rprintf("%s", "Coefficients:\n");
	int z_coef_first = z_coef_ind[0];
	int z_coef_last = z_coef_ind[z_coef_ind.size() - 1];
	NumericMatrix z_coef_results = results(Range(z_coef_first, z_coef_last), _);
	rownames(z_coef_results) = results_rownames[z_coef_ind];
	colnames(z_coef_results) = results_colnames;
	print(as_table(cbind(z_coef_results, stars[z_coef_ind])));

	Rprintf("%s", "---\n");

	Rprintf("%s", "Distribution parameters:\n");
	int distr_first = 0;
	int distr_last = df - z_coef_ind.size() - 1;
	NumericMatrix distr_results = results(Range(distr_first, distr_last), _);
	StringVector distr_rownames = results_rownames[Range(distr_first, distr_last)];
	rownames(distr_results) = distr_rownames;
	colnames(distr_results) = results_colnames;
	print(as_table(cbind(distr_results, stars[Range(distr_first, distr_last)])));

	Rprintf("%s", "---\n");

	Rprintf("%s", "Fixed Coefficients:\n");
	if (is_z_constant_fixed)
	{
		std::string new_str = "(Intercept) = " + std::to_string(z_constant_fixed) + "\n";
		Rprintf("%s", new_str.c_str());
	}

	if (is_z_coef_first_fixed)
	{
		String new_str_names = data_names(1);
		std::string new_str_names_str = new_str_names;
		std::string new_str = new_str_names_str + " = 1" + "\n";
		Rprintf("%s", new_str.c_str());
	}

	Rprintf("%s", "---\n");

	Rprintf("%s", "Fixed Distribution Parameters:\n");
	Rprintf("%s", "a_0 = 1\n");
	if (is_z_mean_fixed)
	{
		std::string new_str = "mean = " + std::to_string(mean) + "\n";
		Rprintf("%s", new_str.c_str());
	}

	if (is_z_sd_fixed)
	{
		std::string new_str = "sd = " + std::to_string(sd) + "\n";
		Rprintf("%s", new_str.c_str());
	}

	Rprintf("%s", "---\n");
	Rprintf("%s", "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n");

	Rprintf("%s", "--------------------------------------------------------------\n");
}

//' Plot hpaBinary random errors approximated density
//' @param x Object of class "hpaBinary"
//' @export	
// [[Rcpp::export]]
void plot_hpaBinary(List x) {

	List model = x;

	//Load additional environments
	Rcpp::Environment base_env("package:base");
	Rcpp::Function round_R = base_env["round"];

	Rcpp::Environment graphics_env("package:graphics");
	Rcpp::Function plot_R = graphics_env["plot"];

	//Load data from the model

	NumericVector pol_coefficients = model["pol_coefficients"];
	NumericVector pol_degrees = model["pol_degrees"];

	NumericVector mean = model["mean"];
	NumericVector sd = model["sd"];

	//Adjust precision

	double errors_exp = model["errors_exp"];
	double errors_var = model["errors_var"];

	double plot_min = errors_exp - 3 * sqrt(errors_var);
	double plot_max = errors_exp + 3 * sqrt(errors_var);

	int n = 10000;

	double precise = (plot_max - plot_min) / n;

	NumericMatrix x_matr = NumericMatrix(n, 1);
	x_matr(0, 0) = plot_min;
	
	for (int i = 1; i < n; i++)
	{
		x_matr(i, 0) = x_matr(i - 1, 0) + precise;
	}

	NumericVector x_vec = x_matr(_, 0);

	//Calculate densities

	NumericVector den = dhpa(x_matr,
		pol_coefficients, pol_degrees,
		LogicalVector::create(false),
		LogicalVector::create(false),
		mean, sd);

	double den_min = min(den);
	double den_max = max(den);

	//Build the plot

	plot_R(Rcpp::_["x"] = x_vec, Rcpp::_["y"] = den,
		Rcpp::_["xlim"] = NumericVector::create(plot_min, plot_max),
		Rcpp::_["xaxp"] = NumericVector::create(plot_min, plot_max, 5),
		Rcpp::_["yaxp"] = NumericVector::create(den_min, den_max, 5),
		Rcpp::_["type"] = "l",
		Rcpp::_["lwd"] = 3,
		Rcpp::_["main"] = "Random Errors Density Approximation Plot",
		Rcpp::_["xlab"] = "value", 
		Rcpp::_["ylab"] = "density");
}

//' Calculates AIC for "hpaBinary" object
//' @description This function calculates AIC for "hpaBinary" object
//' @param object Object of class "hpaBinary"
//' @template AIC_Template
//' @export	
// [[Rcpp::export]]
double AIC_hpaBinary(List object, double k = 2)
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

//' Calculates log-likelihood for "hpaBinary" object
//' @description This function calculates log-likelihood for "hpaBinary" object
//' @param object Object of class "hpaBinary"
//' @export	
// [[Rcpp::export]]
double logLik_hpaBinary(List object) 
{

	double lnL = object["log-likelihood"];

	return(lnL);
}
