#include "hpaMain.h"
#include "hpaML.h"
#include "polynomialIndex.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends(RcppArmadillo)]]

//' Semi-nonparametric maximum likelihood estimation
//' @description This function performs semi-nonparametric maximum likelihood estimation
//' via hermite polynomial densities approximation.
//' @template x_ML_Template
//' @template pol_degrees_Template
//' @template tr_left_Template
//' @template tr_right_Template
//' @template given_ind_Template
//' @template omit_ind_Template
//' @template x0_ML_Template
//' @template hpa_likelihood_details_Template
//' @template GN_details_Template
//' @template first_coef_Template
//' @template parametric_paradigm_Template
//' @template optim_details_Template
//' @return This function returns an object of class "hpaML".\cr \cr
//' An object of class "hpaML" is a list containing the following components:
//' \itemize{
//' \item \code{optim} - \code{\link[stats]{optim}} function output.
//' \item \code{x1} - numeric vector of distribution parameters estimates.
//' \item \code{mean} - density function mean vector estimate.
//' \item \code{sd} - density function sd vector estimate.
//' \item \code{pol_coefficients} - polynomial coefficients estimates.
//' \item \code{tr_left }- the same as \code{tr_left} input parameter.
//' \item \code{tr_right} - the same as \code{tr_right} input parameter.
//' \item \code{omit_ind }- the same as \code{omit_ind} input parameter.
//' \item \code{given_ind} - the same as \code{given_ind} input parameter.
//' \item \code{cov_matrix} - estimated parameters covariance matrix estimate.
//' \item \code{results} - numeric matrix representing estimation results.
//' \item \code{log-likelihood} - value of Log-Likelihood function.
//' \item \code{AIC} - AIC value.
//' \item \code{data} - the same as \code{x} input parameter but without \code{NA} observations.
//' \item \code{n_obs} - number of observations.}
//' @seealso \link[hpa]{summary.hpaML}, \link[hpa]{predict.hpaML}, \link[hpa]{AIC.hpaML}, \link[hpa]{logLik.hpaML}
//' @template hpaML_examples_Template
//' @export
// [[Rcpp::export]]
List hpaML(NumericMatrix x,
	NumericVector pol_degrees = NumericVector(0),
	NumericVector tr_left = NumericVector(0),
	NumericVector tr_right = NumericVector(0),
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector x0 = NumericVector(0))
{

	//Load additional environments

	//stats environment
	Rcpp::Environment stats_env("package:stats");
	Rcpp::Function optim = stats_env["optim"];
	Rcpp::Function na_omit_R = stats_env["na.omit"];

	//base environment
	Rcpp::Environment base_env("package:base");
	Rcpp::Function solve = base_env["solve"];
	Rcpp::Function c_R = base_env["c"];

	//Remove na values from inout vector
	x = na_omit_R(x);

	//Initialize additional index for some loops

	int k = 0;
	
	//Initialize polynomial structure related values

	int pol_degrees_n = pol_degrees.size();

	int pol_coefficients_n = 1;

	for (int i = 0; i < pol_degrees_n; i++)
	{
		pol_coefficients_n *= (pol_degrees[i] + 1);
	}

	pol_coefficients_n -= 1; //because a(0...0)=1

	//Initialize conditions and marginals

	if (given_ind.size() == 0)
	{
		given_ind = LogicalVector(pol_degrees_n);
	}

	if (omit_ind.size() == 0)
	{
		omit_ind = LogicalVector(pol_degrees_n);
	}

	//Determine whether initial values have been manually provided
	bool x0_given = true;

	if (x0.size() == 0)
	{
		x0_given = false;
		x0 = NumericVector(pol_coefficients_n + 2 * pol_degrees_n);
	}

	//Assign indexes

		//For polynomial coefficients
	
	NumericVector pol_coefficients_ind(pol_coefficients_n);

	for (int i = 0; i < pol_coefficients_n; i++)
	{
		pol_coefficients_ind[i] = i;
	}

		//For mean vector

	NumericVector mean_ind(pol_degrees_n);

	for (int i = pol_coefficients_n; i < (pol_coefficients_n + pol_degrees_n); i++)
	{
		mean_ind[k] = i;
		if (!x0_given)
		{
			x0[i] = mean(x(_, k));
		}
		k++;
	}
	
		//For sd vector

	k = 0;

	NumericVector sd_ind(pol_degrees_n);

	for (int i = (pol_coefficients_n + pol_degrees_n); i < (pol_coefficients_n + 2 * pol_degrees_n); i++)
	{
		sd_ind[k] = i;
		if (!x0_given)
		{
			x0[i] = sd(x(_, k));
		}
		k++;
	}

	//Deal with truncation

	NumericMatrix tr_left_mat(1,pol_degrees_n);
	NumericMatrix tr_right_mat(1,pol_degrees_n);

	if ((tr_left.size() > 0) | (tr_right.size() > 0))
	{
		if (tr_left.size() == 0)
		{
			tr_left = NumericVector(pol_degrees_n, R_NegInf);
		}

		if (tr_right.size() == 0)
		{
			tr_right = NumericVector(pol_degrees_n, R_PosInf);
		}

		tr_left_mat(0, _) = tr_left;
		tr_right_mat(0, _) = tr_right;
	} else {
		std::fill(tr_left_mat.begin(), tr_left_mat.end(), NA_REAL);
		std::fill(tr_right_mat.begin(), tr_right_mat.end(), NA_REAL);
	}
	
	//Apply optimization routine

	List PGN_control = List::create(Named("maxit") = 100000);

	List optim_results = optim(
		Rcpp::_["par"] = x0,
		Rcpp::_["fn"] = Rcpp::InternalFunction(&hpaLnLOptim),
		Rcpp::_["control"] = PGN_control,
		Rcpp::_["method"] = "BFGS",
		Rcpp::_["hessian"] = true,
		Rcpp::_["x"] = x,
		Rcpp::_["pol_coefficients_ind"] = pol_coefficients_ind,
		Rcpp::_["pol_degrees"] = pol_degrees,
		Rcpp::_["given_ind"] = given_ind,
		Rcpp::_["omit_ind"] = omit_ind,
		Rcpp::_["mean_ind"] = mean_ind,
		Rcpp::_["sd_ind"] = sd_ind,
		Rcpp::_["is_minus"] = true, //true because of minimization
		Rcpp::_["tr_left"] = tr_left_mat,
		Rcpp::_["tr_right"] = tr_right_mat); 

	//Extract optimization results and assign them to the variables
	//representing estimated parameters

	NumericVector x1 = optim_results[0];
	int x1_n = x1.size();

	double lnL = optim_results["value"];
	lnL = lnL * (-1);

	NumericVector mean = x1[mean_ind];
	NumericVector sd = x1[sd_ind];

	NumericVector pol_coefficients = x1[pol_coefficients_ind];
	pol_coefficients.push_front(1);

	NumericMatrix cov_mat = NumericMatrix(x1_n, x1_n);

	try
	{
		cov_mat = solve(optim_results["hessian"]);
	}
	catch (std::exception &ex)
	{
		warning("Can't estimate covariance matrix because of hessian singulatiry. All covariance matrix elements will be coerced to zero. Hence p-values and significance levels are not interpretable.");
	}
	//Prepare beautifull results output

	NumericMatrix results(x1_n, 3);

	StringVector results_cols = StringVector::create("Estimate", "Std. Error", "P(>|z|)");
	StringVector results_rows(x1_n);

		//Get vector index matrix for polynomial coefficients

	NumericMatrix pol_ind = polynomialIndex(pol_degrees);

	//Assign results matrix columns with values

		//For polynomial coefficients

	for (int i = 1; i < (pol_coefficients_n + 1); i++)
	{
		results_rows[(i - 1)] = "a";
		for (int j = 0; j < pol_degrees_n; j++)
		{
			//Convert double to string
			int my_int = pol_ind(j, i);
			results_rows[(i - 1)] = as<std::string>(results_rows[(i - 1)]) + "_" + std::to_string(my_int);
		}
		results((i - 1), 0) = pol_coefficients[i];
		results((i - 1), 1) = sqrt(cov_mat((i - 1), (i - 1)));
		double t_stat = results((i - 1), 0) / results((i - 1), 1);
		NumericVector F_t_stat = pnorm(NumericVector::create(t_stat));
		results((i - 1), 2) = 2 * std::min(F_t_stat[0], 1 - F_t_stat[0]);
	}

		//For mean

	for (int i = 0; i < pol_degrees_n; i++)
	{
		if (pol_degrees_n == 1)
		{
			results_rows[mean_ind[i]] = "mean";
		} else {
			results_rows[mean_ind[i]] = "mean_" + std::to_string(i + 1);
		}
		results(mean_ind[i], 0) = x1[mean_ind[i]];
		results(mean_ind[i], 1) = sqrt(cov_mat((mean_ind[i]), (mean_ind[i])));
		double z_stat = results(mean_ind[i], 0) / results(mean_ind[i], 1);
		NumericVector F_z_stat = pnorm(NumericVector::create(z_stat));
		results(mean_ind[i], 2) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}

		//For sd

	for (int i = 0; i < pol_degrees_n; i++)
	{
		if (pol_degrees_n == 1)
		{
			results_rows[sd_ind[i]] = "sd";
		} else {
			results_rows[sd_ind[i]] = "sd_" + std::to_string(i + 1);
		}
		results(sd_ind[i], 0) = x1[sd_ind[i]];
		results(sd_ind[i], 1) = sqrt(cov_mat((sd_ind[i]), (sd_ind[i])));
		double z_stat = results(sd_ind[i], 0) / results(sd_ind[i], 1);
		NumericVector F_z_stat = pnorm(NumericVector::create(z_stat));
		results(sd_ind[i], 2) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}

	//Assign names to rows and some vectors

	rownames(results) = results_rows;
	colnames(results) = results_cols;

	x1.names() = results_rows;

	rownames(cov_mat) = results_rows;
	colnames(cov_mat) = results_rows;

	pol_coefficients.names() = c_R("a_0", results_rows[pol_coefficients_ind]);

	//Collect the results
	List return_result = List::create(
		Named("optim") = optim_results, 
		Named("x1") = x1,
		Named("mean") = mean, 
		Named("sd") = sd,
		Named("pol_coefficients") = pol_coefficients, 
		Named("pol_degrees") = pol_degrees,
		Named("tr_left") = tr_left_mat, 
		Named("tr_right") = tr_right_mat,
		Named("omit_ind") = omit_ind, 
		Named("given_ind") = given_ind,
		Named("cov_matrix") = cov_mat,
		Named("results") = results, 
		Named("log-likelihood") = lnL, 
		Named("AIC") = 2 * (x1_n - lnL),
		Named("data") = x,
		Named("n_obs") = x.nrow());

	//Assign the class to the output list
	return_result.attr("class") = "hpaML";

	return(return_result);
}

//Perform log-likelihood function estimation for Phillips-Gallant-Nychka distribution at point
double hpaLnLOptim(NumericVector x0,
	NumericMatrix x = NumericMatrix(1, 1),
	NumericVector pol_coefficients_ind = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector mean_ind = NumericVector(0),
	NumericVector sd_ind = NumericVector(0),
	bool is_minus = false,
	NumericMatrix tr_left = NumericMatrix(1, 1),
	NumericMatrix tr_right = NumericMatrix(1, 1))
{

	NumericVector mean = x0[mean_ind];

	NumericVector sd = x0[sd_ind];

	NumericVector pol_coefficients = x0[pol_coefficients_ind];
	pol_coefficients.push_front(1);

	if (!(R_IsNA(tr_left(0, 0))) & !(R_IsNA(tr_right(0, 0))))
	{

		return((1 - 2 * is_minus) * sum(log(dtrhpa(x,
			tr_left, tr_right,
			pol_coefficients, pol_degrees,
			given_ind, omit_ind,
			mean, sd))));
	}

	return((1 - 2 * is_minus) * sum(log(dhpa(x,
		pol_coefficients, pol_degrees,
		given_ind, omit_ind,
		mean, sd))));
}

//' Predict method for hpaML
//' @param object Object of class "hpaML"
//' @template newdata_Template
//' @return This function returns predictions based on \code{\link[hpa]{hpaML}} estimation results.
//' @export
// [[Rcpp::export]]
NumericVector predict_hpaML(List object, NumericMatrix newdata = NumericMatrix(1,1))
{

	List model = object;

	//Load additional environments

	//stats environment
	Rcpp::Environment stats_env("package:stats");
	Rcpp::Function na_omit_R = stats_env["na.omit"];

	//get distribution parameters

	NumericVector pol_coefficients = model["pol_coefficients"];
	NumericVector pol_degrees = model["pol_degrees"];

	NumericVector mean = model["mean"];
	NumericVector sd = model["sd"];

	NumericMatrix tr_left = model["tr_left"];
	NumericMatrix tr_right = model["tr_right"];

	LogicalVector omit_ind = model["omit_ind"];
	LogicalVector given_ind = model["given_ind"];

	NumericMatrix x = model["data"];

	//get data
	if ((newdata.ncol() == 1) & (newdata.nrow() == 1))
	{
		newdata = x;
	} else {
		newdata = na_omit_R(newdata);
	}

	//estimate
	if (!(R_IsNA(tr_left(0, 0))) & !(R_IsNA(tr_right(0, 0))))
	{
		return(dtrhpa(newdata,
			tr_left, tr_right,
			pol_coefficients, pol_degrees,
			given_ind, omit_ind,
			mean, sd));
	}

	return(dhpa(newdata,
		pol_coefficients, pol_degrees,
		given_ind, omit_ind,
		mean, sd));
}

//' Summarizing hpaML Fits
//' @param object Object of class "hpaML"
//' @return This function returns the same list as \code{\link[hpa]{hpaML}} function changing it's class to "summary.hpaML".
//' @export
// [[Rcpp::export]]
List summary_hpaML(List object)
{

	List return_result = clone(object); //in order to preserve model class

	return_result.attr("class") = "summary.hpaML";

	return(return_result);
}

//' Summary for hpaML output
//' @param x Object of class "hpaML"
//' @export
// [[Rcpp::export]]
void print_summary_hpaML(List x)
{

	List model = x;

	//Load additional environments
	Rcpp::Environment base_env("package:base");
	Rcpp::Function as_table = base_env["as.table"];
	Rcpp::Function cbind = base_env["cbind"];
	Rcpp::Function round_R = base_env["round"];

	NumericMatrix results = model["results"];
	results = round_R(Rcpp::_["x"] = results, Rcpp::_["digits"] = 5);

	NumericVector x1 = model["x1"];

	NumericVector p_values = results(_, 2);

	StringVector stars = starVector(p_values);

	double lnL = model["log-likelihood"];
	double AIC = model["AIC"];
	int n_obs = model["n_obs"];
	int df = x1.size();
	std::string lnL_string = "Log-Likelihood: " + std::to_string(lnL) + "\n";
	std::string AIC_string = "AIC: " + std::to_string(AIC) + "\n";
	std::string n_obs_string = "Observations: " + std::to_string(n_obs) + "\n";
	std::string df_string = std::to_string(df) + " free parameters (df = " + std::to_string(n_obs - df) + ")" + "\n";

	Rprintf("%s", "--------------------------------------------------------------\n");

	Rprintf("%s", "Semi-nonparametric maximum likelihood estimation\n");

	std::printf("%s", lnL_string.c_str());
	std::printf("%s", AIC_string.c_str());
	std::printf("%s", n_obs_string.c_str());
	std::printf("%s", df_string.c_str());

	Rprintf("%s", "Distribution parameters:\n");
	print(as_table(cbind(results, stars)));

	Rprintf("%s", "---\n");
	Rprintf("%s", "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n");

	Rprintf("%s", "--------------------------------------------------------------\n");
}

// Create star vector
StringVector starVector(NumericVector p_values)
{
	int n = p_values.size();

	StringVector stars(n);

	for (int i = 0; i < n; i++)
	{
		if (!NumericVector::is_na(p_values[i]))
		{
			if (p_values[i] <= 0.001)
			{
				stars[i] = "***";
			} else
			{
				if ((0.001 < p_values[i]) & (p_values[i] <= 0.01))
				{
					stars[i] = "**";
				} else
				{
					if ((0.01 < p_values[i]) & (p_values[i] <= 0.05))
					{
						stars[i] = "*";
					} else
					{
						if ((0.05 < p_values[i]) & (p_values[i] <= 0.1))
						{
							stars[i] = ".";
						} else
						{
							stars[i] = " ";
						}
					}
				}
			}
		}
		else
		{
			stars[i] = " ";
		}
	}

	return(stars);
}

//' Calculates AIC for "hpaML" object
//' @description This function calculates AIC for "hpaML" object
//' @param object Object of class "hpaML"
//' @template AIC_template
//' @export
// [[Rcpp::export]]
double AIC_hpaML(List object, double k = 2)
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

//' Calculates log-likelihood for "hpaML" object
//' @description This function calculates log-likelihood for "hpaML" object
//' @param object Object of class "hpaML"
//' @export
// [[Rcpp::export]]
double logLik_hpaML(List object)
{

	double lnL = object["log-likelihood"];

	return(lnL);
}
