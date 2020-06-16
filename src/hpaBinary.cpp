#include "hpaMain.h"
#include "hpaML.h"
#include "hpaBinary.h"
#include "polynomialIndex.h"
#include "ParallelFunctions.h"
#include <RcppParallel.h>

#include <RcppArmadillo.h>
using namespace RcppArmadillo;
using namespace RcppParallel;

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
//' @template cov_type_Template
//' @template boot_iter_Template
//' @template is_parallel_Template
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
//' \item \code{marginal_effects} - marginal effects matrix where columns are variables and rows are observations.
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
	NumericVector x0 = NumericVector(0),
	String cov_type = "sandwich",
	int boot_iter = 100,
	bool is_parallel = false) 
{
  
	// Load additional environments
	
	  // stats environment
	Rcpp::Environment stats_env("package:stats");
	Rcpp::Function optim = stats_env["optim"];
	Rcpp::Function model_frame = stats_env["model.frame"];
	Rcpp::Function na_omit_R = stats_env["na.omit"];
	Rcpp::Function glm = stats_env["glm"];
	Rcpp::Function binomial = stats_env["binomial"];
	Rcpp::Function na_pass = stats_env["na.pass"];
	Rcpp::Function cov_R = stats_env["cov"];
	
	  // base environment
	Rcpp::Environment base_env("package:base");
	Rcpp::Function class_R = base_env["class"];
	Rcpp::Function c_R = base_env["c"];
	Rcpp::Function diag_R = base_env["diag"];
	
	// Initialize polynomial structure related values

	int pol_coefficients_n = K;

	// Initialize bool values related to fixed parameters

	bool is_z_mean_fixed = !R_IsNA(z_mean_fixed);

	bool is_z_sd_fixed = !R_IsNA(z_sd_fixed);

	bool is_z_constant_fixed = !R_IsNA(z_constant_fixed);

	// Working with Data

		// Extract dataframe from formula

	DataFrame z_df = model_frame(Rcpp::_["formula"] = formula, 
                               Rcpp::_["data"] = data,
                               Rcpp::_["na.action"] = na_pass);

	z_df = na_omit_R(z_df);

	int z_df_n = z_df.size();

		// Extract binary dependend variable
	 
	NumericVector z = z_df[0]; // it is reference

	int n_obs = z.size();

		// Extract independend variables

	NumericMatrix z_d(n_obs, (z_df_n - 1) + !is_z_constant_fixed); // -1 because of dependent variable

	int z_d_col = z_d.ncol();

	    // the constant located in last column of regressors matrix
	if (!is_z_constant_fixed)
	{
		z_d(_, z_d_col -  1) = (NumericVector(n_obs) + 1); // add constant

	}

	for (int i = 0; i < (z_d_col - !is_z_constant_fixed); i++)
	{
		z_d(_, i) = as<NumericVector>(z_df[i + 1]); // +1 because of dependent variable
	}

	// Sequential estimation
	if (is_sequence)
	{
		// for K=0
		List results(K + 1);
		results[0] = hpaBinary(formula, data, 
                         0, z_mean_fixed, z_sd_fixed, z_constant_fixed, 
                         is_z_coef_first_fixed, true, false, NumericVector(0), 
                         "sandwich", 100, is_parallel);
		List results_current = results[0];
		NumericVector x1 = results_current["x1"];
		int x0_n = x1.size() + 1; // add one more alpha parameter for the next estimation
		NumericVector x0 = NumericVector(x0_n);
		for (int i = 1; i < x0_n; i++)
		{
			x0[i] = x1[i - 1];
		}
		// for other K
		for (int i = 1; i <= K; i++)
		{
			if (is_z_sd_fixed)
			{
				z_sd_fixed = results_current["sd"];
			}
			results[i] = hpaBinary(formula, data, i, 
                          z_mean_fixed, z_sd_fixed, z_constant_fixed, 
                          is_z_coef_first_fixed, 
                          false, false, x0, "sandwich", 100, is_parallel);
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

	// Create initial values vector
	bool x0_given = true;

	if (x0.size() == 0)
	{
		x0_given = false;
		// x0 dimensions are equal to the estimated (nonfixed) parameters number
		// note thet constant is already in z_d_col or not
		x0 = NumericVector(pol_coefficients_n +
							         !is_z_mean_fixed + !is_z_sd_fixed +
							         z_d_col - is_z_coef_first_fixed);
	}

	// Assign indexes

		// Initialize additional index and upper value for some loops

	int k = 0; // to account for fixed parameters

	int lower_ind = 0;

	int upper_ind = pol_coefficients_n;

		// for polynomial coefficients

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
		// for mean vector

	int z_mean_ind = 0;

	if (!is_z_mean_fixed)
	{
		z_mean_ind = pol_coefficients_n + k;
		k++;
	}

		// for sd vector
	int z_sd_ind = 0;

	if (!is_z_sd_fixed)
	{
		z_sd_ind = pol_coefficients_n + k;
		k++;
	}
	
		// for z coefficients
		// note that z_d_col may contain or not the constant term
	int n_coef = z_d_col - is_z_coef_first_fixed; // number of estimated coefficients

	NumericVector z_coef_ind(n_coef);

	lower_ind = pol_coefficients_n + k;

	upper_ind = pol_coefficients_n + z_coef_ind.size() + k;

	for (int i = lower_ind; i < upper_ind; i++)
	{
		z_coef_ind[i - lower_ind] = i;
	}

	// Convert to arma

	arma::vec z_arma = as<arma::vec>(z);

	arma::mat z_d_arma = as<arma::mat>(z_d);

	// Divide into 0 and 1 samples

	arma::vec z_1 = as<arma::vec>(z[z == 1]);

	arma::mat z_d_1 = (as<arma::mat>(z_d)).rows(arma::find(z_arma == 1));

	arma::vec z_0 = as<arma::vec>(z[z == 0]);

	arma::mat z_d_0 = (as<arma::mat>(z_d)).rows(arma::find(z_arma == 0));

	// set initial sd value to 1

	if (!is_z_sd_fixed & !x0_given)
	{
		x0[z_sd_ind] = 1;
	}

	// Estimate intial values via probit model using glm function

	if (is_x0_probit & !x0_given)
	{
	  // calculate probit model
		List model_probit = glm(Rcpp::_["formula"] = formula, Rcpp::_["data"] = data,
			Rcpp::_["family"] = binomial(Rcpp::_["link"] = "probit"));
	  
	  // extract probit model coefficients estimates
		NumericVector glm_coef = model_probit["coefficients"];
		
		// coefficient for the first regressor which under some
		// input parameters should be used for adjust purposes
		double coef_adjust = std::abs(glm_coef[1]);
		
		if (is_z_coef_first_fixed)
		{
		  // adjust all coefficients
		  glm_coef = glm_coef / coef_adjust;
		  
			// addjust sd to first coefficient value
			if (is_z_sd_fixed)
			{
				z_sd_fixed /= coef_adjust;
			} else {
				x0[z_sd_ind] = x0[z_sd_ind] / coef_adjust;
			}
		} else {
		  // adjust coefficients to sd parameter
			if (is_z_sd_fixed)
			{
				glm_coef = glm_coef / z_sd_fixed;
			}
		}
		if (!is_z_constant_fixed)
		{
			x0[z_coef_ind[n_coef - 1]] = glm_coef[0]; // already adjusted because glm_coef = glm_coef / coef_adjust 
		} else {
			if (!is_z_mean_fixed)
			{
				x0[z_mean_ind] = glm_coef[0]; // already adjusted because glm_coef = glm_coef / coef_adjust 
			} else {
				z_mean_fixed = glm_coef[0];
			}
		}
		for (int i = 0; i < (n_coef - !is_z_constant_fixed); i++)
		{
			x0[z_coef_ind[i]] = glm_coef[i + 1 + is_z_coef_first_fixed]; // + 1 to omitt constant assigned previously
		}
	}

	// Create list for some variables because unfortunatelly optim function has limitation for the
	// parameters number (parameters itself not estimated)

	List is_List = List::create(Named("is_z_coef_first_fixed") = is_z_coef_first_fixed, 
		Named("is_z_mean_fixed") = is_z_mean_fixed,
		Named("is_z_sd_fixed") = is_z_sd_fixed,
		Named("is_z_constant_fixed") = is_z_constant_fixed
	);

	// Apply optimization routine
	
	  // Set optim control parameters

	List PGN_control = List::create(
	     Named("maxit") = 10000000, 
       Named("fnscale") = -1.0,
       Named("abstol") = std::sqrt(std::numeric_limits<double>::epsilon()),
       Named("reltol") = std::sqrt(std::numeric_limits<double>::epsilon()));
	
	    // Perform the optimization

	List optim_results = optim(
		Rcpp::_["par"] = x0,
		Rcpp::_["fn"] = Rcpp::InternalFunction(&hpaBinaryLnLOptim),
		Rcpp::_["gr"] = Rcpp::InternalFunction(&hpaBinaryLnLOptim_grad),
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
	  Rcpp::_["is_parallel"] = is_parallel);

	// Extract coefficients and function value

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

	// Get covariance matrix estimate of "cov_type" type
	
	NumericMatrix cov_mat;
	
	arma::mat H_part;
	
	arma::mat J_part;
	
	NumericMatrix my_hessian;
	
	NumericMatrix my_jacobian;
	
	// Estimate jacobian for the inner part
	
	if ((cov_type == "gop") | (cov_type == "sandwich"))
	{
	  NumericMatrix my_jacobian_tmp = hpaBinaryLnLOptim_grad_ind(x1, is_List,
                                             z_1, z_0,
                                             z_d_1, z_d_0,
                                             K,
                                             z_mean_fixed, z_sd_fixed,
                                             z_constant_fixed, pol_coefficients_ind,
                                             z_mean_ind, z_sd_ind,
                                             z_coef_ind,
                                             is_parallel);
	  
	  my_jacobian = my_jacobian_tmp;
	  
	  J_part = as<arma::mat>(my_jacobian);
	}
	
	if ((cov_type == "hessian") | (cov_type == "sandwich"))
	{
	  NumericMatrix my_hessian_tmp = optim_results["hessian"];
	  my_hessian = my_hessian_tmp;
	  
	  H_part = as<arma::mat>(my_hessian).i();
	}
	
	if ((cov_type == "sandwich"))
	{
	  cov_mat = wrap(H_part * (J_part.t() * J_part) * H_part);
	}
	
	if (cov_type == "gop")
	{
	  cov_mat = wrap((J_part.t() * J_part).i());
	}
	
	if (cov_type == "hessian")
	{
	  cov_mat = wrap(-H_part);
	}
	
	// Apply bootstrap
	
	// store parameters for each iteration
	
	NumericMatrix boot_parameters = NumericMatrix(boot_iter, x1_n);
	
	// store standard deviation
	
	NumericVector sd_dev = NumericVector(x1_n);
	
	// temporal index matrix for each iteration
	
	NumericVector sample_ind = NumericVector(n_obs);
	
	// list to store bootstrap results
	
	List boot_List;
	
	if (cov_type == "bootstrap")
	{
	  for(int i = 0; i < boot_iter; i++)
	  {
	    // Generate sample with replacement
	    
	    NumericVector sample_ind = floor(runif(n_obs, 0, n_obs));
	    
	    NumericVector z_boot = NumericVector(n_obs);
	    
	    NumericMatrix z_d_boot = NumericMatrix(z_d.nrow(), z_d.ncol());
	    
	    
	    for (int j = 0; j < n_obs; j++)
	    {
	      z_boot[j] = z[sample_ind[j]];
	      z_d_boot(j, _) = z_d(sample_ind[j], _);
	    }
	    
	    // Convert to arma
	    
	    arma::vec z_arma_boot = as<arma::vec>(z_boot);
	    
	    arma::mat z_d_arma_boot = as<arma::mat>(z_d_boot);
	    
	    // Divide into 0 and 1 samples
	    
	    arma::vec z_1_boot = as<arma::vec>(z_boot[z_boot == 1]);
	    
	    arma::mat z_d_1_boot = (as<arma::mat>(z_d_boot)).rows(arma::find(z_arma_boot == 1));
	    
	    arma::vec z_0_boot = as<arma::vec>(z_boot[z_boot == 0]);
	    
	    arma::mat z_d_0_boot = (as<arma::mat>(z_d_boot)).rows(arma::find(z_arma_boot == 0));
	    
	    // Perform estimaton
	    
	    List boot_results = optim(
	      Rcpp::_["par"] = x1,
	      Rcpp::_["fn"] = Rcpp::InternalFunction(&hpaBinaryLnLOptim),
	      Rcpp::_["gr"] = Rcpp::InternalFunction(&hpaBinaryLnLOptim_grad),
	      Rcpp::_["control"] = PGN_control,
	      Rcpp::_["method"] = "BFGS",
	      Rcpp::_["hessian"] = true,
	      Rcpp::_["is_List"] = is_List,
	      Rcpp::_["z_1"] = z_1_boot,
	      Rcpp::_["z_0"] = z_0_boot,
	      Rcpp::_["z_d_1"] = z_d_1_boot,
	      Rcpp::_["z_d_0"] = z_d_0_boot,
	      Rcpp::_["K"] = K,
	      Rcpp::_["z_mean_fixed"] = z_mean_fixed,
	      Rcpp::_["z_sd_fixed"] = z_sd_fixed,
	      Rcpp::_["z_constant_fixed"] = z_constant_fixed,
	      Rcpp::_["pol_coefficients_ind"] = pol_coefficients_ind,
	      Rcpp::_["z_mean_ind"] = z_mean_ind,
	      Rcpp::_["z_sd_ind"] = z_sd_ind,
	      Rcpp::_["z_coef_ind"] = z_coef_ind,
	      Rcpp::_["is_parallel"] = is_parallel);
	    
	    // Store iteration results
	    
	    NumericVector x1_new = boot_results["par"];
	    
	    boot_parameters(i, _) = x1_new;
	  }
	  
	  // Store bootstrap results
	  
	  cov_mat = cov_R(boot_parameters);
	  
	  sd_dev = sqrt(diag_R(cov_mat));
	  
	  boot_List = List::create(
	    Named("estimates") = boot_parameters,
	    Named("cov_mat") = cov_mat,
	    Named("sd") = sd_dev);
	}

	// Prepare beautifull results output

	NumericMatrix results(x1_n, 4);

	StringVector results_cols = StringVector::create("Estimate", "Std. Error", "z value", "P(>|z|)");

	StringVector results_rows(x1_n);

	NumericMatrix pol_ind = polynomialIndex(NumericVector::create(K));

		// for alpha

		for (int i = 0; i < K; i++)
		{
			// Convert double to string
			std::stringstream ss;
			ss << (i + 1); // to start from alpha_1
			std::string str2 = ss.str();
			//
			results_rows[i] = "a_" + str2;
			results(i, 0) = pol_coefficients[i + 1];
			results(i, 1) = sqrt(cov_mat(i, i));
			double z_stat = results(i, 0) / results(i, 1);
			NumericVector F_z_stat = pnorm(NumericVector::create(z_stat));
			results(i, 2) = z_stat;
			results(i, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
		}

		// for mean
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

		// for sd
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
	
		// for z coefficients
	CharacterVector z_df_names = z_df.names();
	String first_coef_name = z_df_names[1];
	z_df_names.erase(0); // remove dependend variable name
	
	// for intercept if need
	if (!is_z_constant_fixed)
	{
		z_df_names.push_back("(Intercept)");
	}
	
	// remove first regressors if it's coefficient is fixed
	if (is_z_coef_first_fixed)
	{
	  z_df_names.erase(0);
	}

	int coef_start = x1_n - n_coef;
	for (int i = coef_start; i < x1_n; i++)
	{
		results_rows(i) = z_df_names[i - coef_start];//+1 because the first is dependend variable
		results(i, 0) = z_coef[i - coef_start];
		results(i, 1) = sqrt(cov_mat(z_coef_ind[i - coef_start], z_coef_ind[i - coef_start]));
		double z_stat = results(i, 0) / results(i, 1);
		NumericVector F_z_stat = pnorm(NumericVector::create(z_stat));
		results(i, 2) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
		results(i, 3) = 2 * std::min(F_z_stat[0], 1 - F_z_stat[0]);
	}

		// for expectation and variance
	NumericVector z_e = ehpa(NumericMatrix(1, 1), pol_coefficients, NumericVector::create(K),
		LogicalVector::create(false), LogicalVector::create(false),
		z_mean, z_sd, 
		NumericVector::create(1), false);

	NumericVector z_e_2 = ehpa(NumericMatrix(1, 1), pol_coefficients, NumericVector::create(K),
		LogicalVector::create(0), LogicalVector::create(0),
		z_mean, z_sd,
		NumericVector::create(2), false);

		// assign names to the output
	rownames(results) = results_rows;

	colnames(results) = results_cols;

	x1.names() = results_rows;

	rownames(cov_mat) = results_rows;
	colnames(cov_mat) = results_rows;

	if (K != 0)
	{
		pol_coefficients.names() = c_R("a_0", results_rows[pol_coefficients_ind]);
	}

	StringVector z_coef_names = results_rows[z_coef_ind];

	if (is_z_coef_first_fixed)
	{
	  z_coef.push_front(1); // add identity coefficient for fixed value
	  z_coef_names.push_front(first_coef_name);
	}
	
	z_coef.names() = z_coef_names;

	// Estimate latent variable and probabilities

		// coefficients for independend variables

	arma::vec z_coef_arma = as<arma::vec>(z_coef);

		// get estimates for z*

	NumericMatrix z_latent = wrap(z_d_arma * z_coef_arma);

	if (is_z_constant_fixed)
	{
		z_latent = z_latent + z_constant_fixed;
	}

	NumericVector z_prob = 1 - phpa(-1.0 * z_latent, pol_coefficients,
                              		NumericVector::create(K),
                              		LogicalVector::create(0), 
                              		LogicalVector::create(0),
                              		z_mean, z_sd,
                              		is_parallel);
	
	// Estimate marginal effects

	NumericVector z_den = dhpa(z_latent, pol_coefficients,
                             NumericVector::create(K),
                             LogicalVector::create(0), 
                             LogicalVector::create(0),
                             z_mean, z_sd,
                             is_parallel);
	
	int n_coef_total = z_coef.size();
	
	NumericMatrix marginal_effects = NumericMatrix(z_latent.nrow(), n_coef_total);
	
	for (int i = 0; i < n_coef_total; i++)
	{
	  NumericVector me_vec = z_den * z_coef[i];
	  marginal_effects(_, i) = me_vec;
	}

	StringVector n_coefames = z_coef.names();
	colnames(marginal_effects) = n_coefames;

	// return results
	
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
		Named("log-likelihood") = optim_value,
		Named("AIC") = 2 * (x1_n - optim_value),
		Named("n_obs") = z_latent.nrow(),
		Named("z_latent") = z_latent,
		Named("z_prob") = z_prob,
		Named("formula") = formula,
		Named("dataframe") = z_df,
		Named("model_Lists") = model_Lists,
		Named("cov_mat") = cov_mat,
		Named("marginal_effects") = marginal_effects);

	return_result.attr("class") = "hpaBinary";

		// estimate probabilities 

	return(return_result);
}

// Perform semi-nonparametric log-likelihood function estimation for binary choice model
List hpaBinaryLnLOptim_List(NumericVector x0,
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
	bool is_parallel = false) {

	// Get values from the is_List
	bool is_z_coef_first_fixed = is_List["is_z_coef_first_fixed"];
	bool is_z_mean_fixed = is_List["is_z_mean_fixed"];
	bool is_z_sd_fixed = is_List["is_z_sd_fixed"];
	bool is_z_constant_fixed = is_List["is_z_constant_fixed"];

	// Assign estimated parameters values to corresponding vectors

		// polynomial coefficients and degrees

	NumericVector pol_coefficients = NumericVector(K);

	if (K != 0)
	{
		pol_coefficients = x0[pol_coefficients_ind];
	} 

	pol_coefficients.push_front(1);

	NumericVector pol_degrees = NumericVector(1);

	pol_degrees[0] = K;

		// mean value

	NumericVector z_mean = NumericVector(1);

	if (is_z_mean_fixed)
	{
		z_mean[0] = z_mean_fixed;
	}
	else {
		z_mean[0] = x0[z_mean_ind];
	}

		// sd value

	NumericVector z_sd = NumericVector(1);

	if (is_z_sd_fixed)
	{
		z_sd[0] = z_sd_fixed;
	}
	else {
		z_sd[0] = x0[z_sd_ind];
	}

		// coefficients for independend variables

	NumericVector z_coef_R = x0[z_coef_ind];

	if (is_z_coef_first_fixed)
	{
		z_coef_R.push_front(1); // add identity coefficient for fixed value
	}

	arma::vec z_coef = as<arma::vec>(z_coef_R);

	// get estimates for z*

	NumericMatrix z_h_1 = wrap(z_d_1 * z_coef);
	
	NumericMatrix z_h_0 = wrap(z_d_0 * z_coef);

	  if (is_z_constant_fixed)
  	{
  		z_h_1 = z_h_1 + z_constant_fixed;
	    z_h_0 = z_h_0 + z_constant_fixed;
  	}

	// Likelihood calculation
	
	NumericVector lnL_z_1;
	
	NumericVector lnL_z_0;

  lnL_z_1 = log(1 - phpa(-1.0 * z_h_1,
  		          pol_coefficients, pol_degrees,
  		          LogicalVector{false}, LogicalVector{false},
  		          z_mean, z_sd, 
  		          is_parallel));

  lnL_z_0 = log(phpa(-1.0 * z_h_0,
  		               pol_coefficients, pol_degrees,
  		               LogicalVector{false}, LogicalVector{false},
  		               z_mean, z_sd,
  		               is_parallel));

	// Initialize list to store calculation results
	
	double aggregate_0 = 0.0;
	double aggregate_1 = 0.0;
	
	List return_List = List::create(Named("individual_1") = NumericVector::create(0.0),
                                  Named("individual_0") = NumericVector::create(0.0),
                                  Named("aggregate_1") = aggregate_1,
                                  Named("aggregate_0") = aggregate_0);
	
	// Store calculation results

	return_List["individual_1"] = lnL_z_1;
	aggregate_1 = sum(lnL_z_1);
	return_List["aggregate_1"] = aggregate_1;

	return_List["individual_0"] = lnL_z_0;
	aggregate_0 = sum(lnL_z_0);
	return_List["aggregate_0"] = aggregate_0;
	
	return(return_List);
}

// Perform semi-nonparametric log-likelihood function estimation for binary choice model
double hpaBinaryLnLOptim(NumericVector x0,
                         List is_List,
                         arma::vec z_1, arma::vec z_0,
                         arma::mat z_d_1, arma::mat z_d_0,
                         int K = 1,
                         double z_mean_fixed = NA_REAL,
                         double z_sd_fixed = NA_REAL,
                         double z_constant_fixed = 0,
                         NumericVector pol_coefficients_ind = NumericVector(0),
                         int z_mean_ind = 1,
                         int z_sd_ind = 2,
                         NumericVector z_coef_ind = NumericVector(0),
                         bool is_parallel = false) 
{ 
  List return_List = hpaBinaryLnLOptim_List(x0, 
                                            is_List,
                                            z_1, z_0,
                                            z_d_1, z_d_0,
                                            K,
                                            z_mean_fixed, z_sd_fixed,
                                            z_constant_fixed,
                                            pol_coefficients_ind,
                                            z_mean_ind, z_sd_ind,
                                            z_coef_ind, 
                                            is_parallel);
  
  double aggregate_0 = return_List["aggregate_0"];
  double aggregate_1 = return_List["aggregate_1"];
  
  double return_aggregate = 0.0;
  
  return_aggregate += aggregate_0;
  return_aggregate += aggregate_1;
  
  return(return_aggregate);
}
 
// Perform semi-nonparametric log-likelihood function estimation for binary choice model
NumericVector hpaBinaryLnLOptim_ind(NumericVector x0,
                          List is_List,
                          arma::vec z_1, arma::vec z_0,
                          arma::mat z_d_1, arma::mat z_d_0,
                          int K = 1,
                          double z_mean_fixed = NA_REAL,
                          double z_sd_fixed = NA_REAL,
                          double z_constant_fixed = 0,
                          NumericVector pol_coefficients_ind = NumericVector(0),
                          int z_mean_ind = 1,
                          int z_sd_ind = 2,
                          NumericVector z_coef_ind = NumericVector(0),
                          bool is_parallel = false) 
  { 
   List return_List = hpaBinaryLnLOptim_List(x0, 
                                             is_List,
                                             z_1, z_0,
                                             z_d_1, z_d_0,
                                             K,
                                             z_mean_fixed, z_sd_fixed,
                                             z_constant_fixed,
                                             pol_coefficients_ind,
                                             z_mean_ind, z_sd_ind,
                                             z_coef_ind, 
                                             is_parallel);
   
   NumericVector individual_0 = return_List["individual_0"];
   NumericVector individual_1 = return_List["individual_1"];
   
   int n_obs_0 = individual_0.size();
   int n_obs_1 = individual_1.size();
   int n_obs = n_obs_0 + n_obs_1;
   
   NumericVector return_individual = NumericVector(n_obs);
  
  return_individual[Range(0, n_obs_1 - 1)] = individual_1;
  return_individual[Range(n_obs_1, n_obs - 1)] = individual_0;

  return(return_individual);
 }                                          


List hpaBinaryLnLOptim_grad_List(NumericVector x0,
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
                                 bool is_parallel = false) {
  
  // Get parameters number
  
  int n_param = x0.size();
  
  // Get estimated regressors number
  
  int n_reg = z_coef_ind.size();
  
  // Get values from the is_List
  
  bool is_z_coef_first_fixed = is_List["is_z_coef_first_fixed"];
  
  bool is_z_mean_fixed = is_List["is_z_mean_fixed"];
  
  bool is_z_sd_fixed = is_List["is_z_sd_fixed"];
  
  bool is_z_constant_fixed = is_List["is_z_constant_fixed"];
  
  // Assign estimated parameters values to corresponding vectors
  
    // polynomial coefficients and degrees
  
  NumericVector pol_coefficients = NumericVector(K);
  
  if (K != 0)
  {
    pol_coefficients = x0[pol_coefficients_ind];
  } 
  
  pol_coefficients.push_front(1);
  
  NumericVector pol_degrees = NumericVector(1);
  
  pol_degrees[0] = K;
  
    // mean value
  
  NumericVector z_mean = NumericVector(1);
  
  if (is_z_mean_fixed)
  {
    z_mean[0] = z_mean_fixed;
  }
  else {
    z_mean[0] = x0[z_mean_ind];
  }
  
    // sd value
  
  NumericVector z_sd = NumericVector(1);
  
  if (is_z_sd_fixed)
  {
    z_sd[0] = z_sd_fixed;
  }
  else {
    z_sd[0] = x0[z_sd_ind];
  }
  
    // coefficients for independend variables
  
  NumericVector z_coef_R = x0[z_coef_ind];
  
  if (is_z_coef_first_fixed)
  {
    z_coef_R.push_front(1); // add identity coefficient for fixed value
  }
  
  arma::vec z_coef = as<arma::vec>(z_coef_R);
  
  // get estimates for z*
  
  NumericMatrix z_h_1 = wrap(z_d_1 * z_coef);
  NumericMatrix z_h_0 = wrap(z_d_0 * z_coef);
  
  if (is_z_constant_fixed)
  {
    z_h_1 = z_h_1 + z_constant_fixed;
    z_h_0 = z_h_0 + z_constant_fixed;
  }
  
  // calculate observations numbers
  
  int n_obs_1 = z_h_1.nrow();
  int n_obs_0 = z_h_0.nrow();
  int n_obs = n_obs_0 + n_obs_1;
  
  // cdf calculation
  
  NumericVector cdf_z_1;
  
  NumericVector cdf_z_0;

  cdf_z_1 = 1 - phpa(-1.0 * z_h_1,
                     pol_coefficients, pol_degrees,
                     LogicalVector{false}, LogicalVector{false},
                     z_mean, z_sd,
                     is_parallel);

  cdf_z_0 = phpa(-1.0 * z_h_0,
                 pol_coefficients, pol_degrees,
                 LogicalVector{false}, LogicalVector{false},
                 z_mean, z_sd,
                 is_parallel);
  
  // pdf calculation
  
  NumericVector pdf_z_1;
  
  NumericVector pdf_z_0;

  pdf_z_1 = dhpa(-1.0 * z_h_1,
                 pol_coefficients, pol_degrees,
                 LogicalVector{false}, LogicalVector{false},
                 z_mean, z_sd,
                 is_parallel);

  pdf_z_0 = dhpa(-1.0 * z_h_0,
                 pol_coefficients, pol_degrees,
                 LogicalVector{false}, LogicalVector{false},
                 z_mean, z_sd,
                 is_parallel);
  
  // Initialize vector to store gradient values
  NumericMatrix my_grad = NumericMatrix(n_obs_0 + n_obs_1, n_param);
  
  // Analytical part of gradient
  
    // for polynomial coefficients
  int pol_coefficients_n = pol_coefficients.size();
  
  NumericMatrix pol_coefficients_grad_1;
  NumericMatrix pol_coefficients_grad_0;
  
  // gradient for ones values respect to polynomial coefficients

      // lower tail negative infinity matrix
    NumericMatrix x_lower_1 = NumericMatrix(n_obs_1, 1);
    std::fill(x_lower_1.begin(), x_lower_1.end(), R_NegInf);
      
      //gradient calculations
    pol_coefficients_grad_1 = ihpaDiff(x_lower_1, -1.0 * z_h_1,
                                       pol_coefficients, pol_degrees,
                                       LogicalVector{false}, LogicalVector{false},
                                       z_mean, z_sd,
                                       "pol_coefficients",
                                       is_parallel);
  
  // gradient for zero values respect to polynomial coefficients

      // lower tail negative infinity matrix
    NumericMatrix x_lower_0 = NumericMatrix(n_obs_0, 1);
    std::fill(x_lower_0.begin(), x_lower_0.end(), R_NegInf);
    
      //gradient calculations
    pol_coefficients_grad_0 = ihpaDiff(x_lower_0, -1.0 * z_h_0,
                                       pol_coefficients, pol_degrees,
                                       LogicalVector{false}, LogicalVector{false},
                                       z_mean, z_sd,
                                       "pol_coefficients",
                                       is_parallel);
    
  for (int i = 0; i < (pol_coefficients_n - 1); i++) // for each parameter
  {
    NumericVector my_grad_tmp = NumericVector(n_obs);

    my_grad_tmp[Range(0, n_obs_1 - 1)] = -1.0 * pol_coefficients_grad_1(_, i + 1) / cdf_z_1;
    my_grad_tmp[Range(n_obs_1, n_obs - 1)] = pol_coefficients_grad_0(_, i + 1) / cdf_z_0;

    my_grad(_, i) = my_grad_tmp;
  }

    // for regression coefficients
    
  NumericVector z_coef_grad_shareble_1;
  NumericVector z_coef_grad_shareble_0;
  
  z_coef_grad_shareble_1 = pdf_z_1 / cdf_z_1;
  z_coef_grad_shareble_0 = pdf_z_0 / cdf_z_0;
  
  NumericMatrix z_d_0_NM = wrap(z_d_0);
  NumericMatrix z_d_1_NM = wrap(z_d_1);

  for (int i = 0; i < n_reg; i++) // for each regressor
  {
    NumericVector my_grad_tmp = NumericVector(n_obs);
    
    my_grad_tmp[Range(0, n_obs_1 - 1)] = z_d_1_NM(_, i + is_z_coef_first_fixed) * z_coef_grad_shareble_1;
    my_grad_tmp[Range(n_obs_1, n_obs - 1)] = -1.0 * z_d_0_NM(_, i + is_z_coef_first_fixed) * z_coef_grad_shareble_0;
    
    my_grad(_, z_coef_ind[i]) = my_grad_tmp;
  }

  // Numeric gradients
  
    // Set differentiation increment value for each parameter
  
  double machinePrecision = std::numeric_limits<double>::epsilon();
  double my_precision = std::sqrt(machinePrecision);
  
  NumericVector eps = abs(x0 * my_precision);
  
    // Control for zero values
  
  eps[eps < (machinePrecision * 100)] = my_precision;
  
    // Estimate the gradient itself
  
  NumericVector x0_eps = clone(x0);
  
  NumericVector f1 = NumericVector(n_obs);
  NumericVector f2 = NumericVector(n_obs);

    // for mean
  
  if (!is_z_mean_fixed)
  {
    
      // Calculate f(x - eps)
    
    x0_eps[z_mean_ind] = x0[z_mean_ind] - eps[z_mean_ind];
    
    f1 = hpaBinaryLnLOptim_ind(x0_eps,
                               is_List,
                               z_1, z_0,
                               z_d_1, z_d_0,
                               K,
                               z_mean_fixed, z_sd_fixed,
                               z_constant_fixed,
                               pol_coefficients_ind,
                               z_mean_ind, z_sd_ind,
                               z_coef_ind,
                               is_parallel);
    
      // Calculate f(x - eps)
    
    x0_eps[z_mean_ind] = x0[z_mean_ind] + eps[z_mean_ind];
      
    f2 = hpaBinaryLnLOptim_ind(x0_eps,
                               is_List,
                               z_1, z_0,
                               z_d_1, z_d_0,
                               K,
                               z_mean_fixed, z_sd_fixed,
                               z_constant_fixed,
                               pol_coefficients_ind,
                               z_mean_ind, z_sd_ind,
                               z_coef_ind,
                               is_parallel);
    
      // Estimate the gradient 
    
    my_grad(_, z_mean_ind) = (f2 - f1) / (2 * eps[z_mean_ind]);
    
      // Set x0_eps value to default
    
    x0_eps[z_mean_ind] = x0[z_mean_ind];
  }
  
  if (!is_z_sd_fixed)
  {
    
    // Calculate f(x - eps)
    
    x0_eps[z_sd_ind] = x0[z_sd_ind] - eps[z_sd_ind];
    
    f1 = hpaBinaryLnLOptim_ind(x0_eps,
                               is_List,
                               z_1, z_0,
                               z_d_1, z_d_0,
                               K,
                               z_mean_fixed, z_sd_fixed,
                               z_constant_fixed,
                               pol_coefficients_ind,
                               z_mean_ind, z_sd_ind,
                               z_coef_ind,
                               is_parallel);
    
    // Calculate f(x - eps)
    
    x0_eps[z_sd_ind] = x0[z_sd_ind] + eps[z_sd_ind];
    
    f2 = hpaBinaryLnLOptim_ind(x0_eps,
                               is_List,
                               z_1, z_0,
                               z_d_1, z_d_0,
                               K,
                               z_mean_fixed, z_sd_fixed,
                               z_constant_fixed,
                               pol_coefficients_ind,
                               z_mean_ind, z_sd_ind,
                               z_coef_ind,
                               is_parallel);
    
    // Estimate the gradient 
    
    my_grad(_, z_sd_ind) = (f2 - f1) / (2 * eps[z_sd_ind]);
    
    // Set x0_eps value to default
    
    x0_eps[z_sd_ind] = x0[z_sd_ind];
  }
  
  List return_List = List::create(Named("aggregate") = colSums(my_grad),
                                  Named("individual") = my_grad);
  
  return(return_List);
}

NumericVector hpaBinaryLnLOptim_grad(NumericVector x0,
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
                                         bool is_parallel = false) 
{
  List return_List = hpaBinaryLnLOptim_grad_List(x0, 
                                                 is_List,
                                                 z_1, z_0,
                                                 z_d_1, z_d_0,
                                                 K,
                                                 z_mean_fixed, z_sd_fixed,
                                                 z_constant_fixed,
                                                 pol_coefficients_ind,
                                                 z_mean_ind, z_sd_ind,
                                                 z_coef_ind, 
                                                 is_parallel);
  
  NumericVector return_aggregate = return_List["aggregate"];

  return(return_aggregate);
}

NumericMatrix hpaBinaryLnLOptim_grad_ind(NumericVector x0,
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
                                 bool is_parallel = false) {
  
  List return_List = hpaBinaryLnLOptim_grad_List(x0, 
                                                 is_List,
                                                 z_1, z_0,
                                                 z_d_1, z_d_0,
                                                 K,
                                                 z_mean_fixed, z_sd_fixed,
                                                 z_constant_fixed,
                                                 pol_coefficients_ind,
                                                 z_mean_ind, z_sd_ind,
                                                 z_coef_ind, 
                                                 is_parallel);
  
  NumericMatrix return_individual = return_List["individual"];
  
  return(return_individual);
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

	// Add additional environments

	Rcpp::Environment stats_env("package:stats");
	Rcpp::Function model_frame = stats_env["model.frame"];
	Rcpp::Function na_omit_R = stats_env["na.omit"];

	Rcpp::Environment base_env("package:base");
	Rcpp::Function as_data_frame = base_env["as.data.frame"];

	// Extract variables from model

		// extract is values

	List model_Lists = model["model_Lists"];

	List is_List = model_Lists["is_List"];

	bool is_z_constant_fixed = is_List["is_z_constant_fixed"];

		// extract fixed values

	List fixed_List = model_Lists["fixed_List"];

	double z_constant_fixed = fixed_List["z_constant_fixed"];

		// extract coefficients

	NumericVector pol_coefficients = model["pol_coefficients"];

	NumericVector z_mean = model["mean"];

	NumericVector z_sd = model["sd"];

	NumericVector z_coef = model["coefficients"];

		// extract polynomial coefficients

	double K = pol_coefficients.size() - 1;

	// Check wheather new dataframe has been supplied

	DataFrame data = newdata;

	if (newdata.size() == 0)
	{
		newdata = as_data_frame(model["dataframe"]);
	}

	// Remove NA values

	data = na_omit_R(newdata);

	// Working with Data

		// Extract dataframe from formula

	Formula formula = model["formula"];

	DataFrame z_df = model_frame(Rcpp::_["formula"] = formula, Rcpp::_["data"] = data);

	int z_df_n = z_df.size();

	// Extract binary dependend variable

	NumericVector z = z_df[0]; // it is reference

	int n_obs = z.size();
	
	// Extract independend variables
	
	NumericMatrix z_d(n_obs, (z_df_n - 1) + !is_z_constant_fixed); // -1 because of dependent variable
	
	int z_d_col = z_d.ncol();
	
	  // the constant located in last column of regressors matrix
	if (!is_z_constant_fixed)
	{
	  z_d(_, z_d_col -  1) = (NumericVector(n_obs) + 1); // add constant
	  
	}
	
	for (int i = 0; i < (z_d_col - !is_z_constant_fixed); i++)
	{
	  z_d(_, i) = as<NumericVector>(z_df[i + 1]); // +1 because of dependent variable
	}

		// Convert to arma

	arma::vec z_arma = as<arma::vec>(z);

	arma::mat z_d_arma = as<arma::mat>(z_d);

	// Estimate latent variable and probabilities

		// coefficients for independend variables

	arma::vec z_coef_arma = as<arma::vec>(z_coef);

	// get estimates for z*

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
		z_mean, z_sd, false);

	return(z_prob);
}

//' Summarizing hpaBinary Fits
//' @param object Object of class "hpaBinary"
//' @return This function returns the same list as \code{\link[hpa]{hpaBinary}} function changing it's class to "summary.hpaBinary".
//' @export
// [[Rcpp::export]]
List summary_hpaBinary(List object)
{

	List return_result = clone(object); // in order to preserve model class

	return_result.attr("class") = "summary.hpaBinary";

	return(return_result);
}

//' Summary for hpaBinary output
//' @param x Object of class "hpaML"
//' @export	
// [[Rcpp::export]]
void print_summary_hpaBinary(List x)
{

	// Extract the model

	List model = x;

	// Load additional environments
	Rcpp::Environment base_env("package:base");
	Rcpp::Function as_table = base_env["as.table"];
	Rcpp::Function cbind = base_env["cbind"];
	Rcpp::Function round_R = base_env["round"];
	Rcpp::Function as_data_frame = base_env["as.data.frame"];
	Rcpp::Function print_R = base_env["print"];
	Rcpp::Function cat_R = base_env["cat"];

	NumericMatrix results = model["results"];
	results = round_R(Rcpp::_["x"] = results, Rcpp::_["digits"] = 5);

	List model_Lists = model["model_Lists"];

	List ind_List = model_Lists["ind_List"];
	NumericVector z_coef_ind = ind_List["z_coef_ind"];

	// Extract is values

	List is_List = model_Lists["is_List"];
	bool is_z_coef_first_fixed = is_List["is_z_coef_first_fixed"];
	bool is_z_mean_fixed = is_List["is_z_mean_fixed"];
	bool is_z_sd_fixed = is_List["is_z_sd_fixed"];
	bool is_z_constant_fixed = is_List["is_z_constant_fixed"];

	// Extract fixed values

	List fixed_List = model_Lists["fixed_List"];
	double z_constant_fixed = fixed_List["z_constant_fixed"];

	// Other stuff

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

	cat_R("--------------------------------------------------------------\n");

	cat_R("Semi-nonparametric binary choice model estimation\n");

	cat_R(lnL_string.c_str());
	cat_R(AIC_string.c_str());
	cat_R(n_obs_string.c_str());
	cat_R(df_string.c_str());
	
	cat_R("---\n");

	cat_R("Coefficients:\n");
	int z_coef_first = z_coef_ind[0];
	int z_coef_last = z_coef_ind[z_coef_ind.size() - 1];
	NumericMatrix z_coef_results = results(Range(z_coef_first, z_coef_last), _);
	rownames(z_coef_results) = results_rownames[z_coef_ind];
	colnames(z_coef_results) = results_colnames;
	print_R(as_table(cbind(z_coef_results, stars[z_coef_ind])));

	cat_R("---\n");

	cat_R("Distribution parameters:\n");
	int distr_first = 0;
	int distr_last = df - z_coef_ind.size() - 1;
	NumericMatrix distr_results = results(Range(distr_first, distr_last), _);
	StringVector distr_rownames = results_rownames[Range(distr_first, distr_last)];
	rownames(distr_results) = distr_rownames;
	colnames(distr_results) = results_colnames;
	print(as_table(cbind(distr_results, stars[Range(distr_first, distr_last)])));

	cat_R("---\n");

	cat_R("Fixed Coefficients:\n");
	if (is_z_constant_fixed)
	{
		std::string new_str = "(Intercept) = " + std::to_string(z_constant_fixed) + "\n";
	  cat_R(new_str.c_str());
	}

	if (is_z_coef_first_fixed)
	{
		String new_str_names = data_names(1);
		std::string new_str_names_str = new_str_names;
		std::string new_str = new_str_names_str + " = 1" + "\n";
		cat_R(new_str.c_str());
	}

	cat_R("---\n");

	cat_R("Fixed Distribution Parameters:\n");
	cat_R("a_0 = 1\n");
	if (is_z_mean_fixed)
	{
		std::string new_str = "mean = " + std::to_string(mean) + "\n";
	  cat_R(new_str.c_str());
	}

	if (is_z_sd_fixed)
	{
		std::string new_str = "sd = " + std::to_string(sd) + "\n";
	  cat_R(new_str.c_str());
	}

	cat_R("---\n");
	cat_R("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n");

	cat_R("--------------------------------------------------------------\n");
}

//' Plot hpaBinary random errors approximated density
//' @param x Object of class "hpaBinary"
//' @export	
// [[Rcpp::export]]
void plot_hpaBinary(List x) {

	List model = x;

	// Load additional environments
	Rcpp::Environment base_env("package:base");
	Rcpp::Function round_R = base_env["round"];
	Rcpp::Function plot_R = base_env["plot"];

	// Load data from the model

	NumericVector pol_coefficients = model["pol_coefficients"];
	NumericVector pol_degrees = model["pol_degrees"];

	NumericVector mean = model["mean"];
	NumericVector sd = model["sd"];

	// Adjust precision

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

	// Calculate densities

	NumericVector den = dhpa(x_matr,
		pol_coefficients, pol_degrees,
		LogicalVector::create(false),
		LogicalVector::create(false),
		mean, sd, false);

	double den_min = min(den);
	double den_max = max(den);

	// Build the plot

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
