#include "hpaMain.h"
#include "hpaML.h"
#include "polynomialIndex.h"
#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace RcppParallel;
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
//' @template cov_type_Template
//' @template boot_iter_Template
//' @template is_parallel_Template
//' @template hpa_likelihood_details_Template
//' @template GN_details_Template
//' @template first_coef_Template
//' @template parametric_paradigm_Template
//' @template optim_details_Template
//' @return This function returns an object of class "hpaML".\cr \cr
//' An object of class "hpaML" is a list containing the following components:
//' \itemize{
//' \item \code{optim} - \code{\link[maxLik]{maxLik}} function output.
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
//' \item \code{n_obs} - number of observations.
//' \item \code{bootstrap} - list where bootstrap estimation results are stored.}
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
	NumericVector x0 = NumericVector(0),
	String cov_type = "sandwich",
	int boot_iter = 100,
	bool is_parallel = false)
{

	// Load additional environments

		// stats environment

	Rcpp::Environment stats_env("package:stats");
	Rcpp::Function na_omit_R = stats_env["na.omit"];
	Rcpp::Function optim = stats_env["optim"];
	Rcpp::Function cov_R = stats_env["cov"];

		// base environment

	Rcpp::Environment base_env("package:base");
	Rcpp::Function c_R = base_env["c"];
	Rcpp::Function diag_R = base_env["diag"];
  
	// Remove NA values from data

	x = na_omit_R(x);
	
	// Get the number of observations

	int n_obs = x.nrow();
	
	// Initialize polynomial structure related values

	int pol_degrees_n = pol_degrees.size();       // random vector dimensionality

	int pol_coefficients_n = 1;                   // number of polynomial coefficients

	for (int i = 0; i < pol_degrees_n; i++)
	{
		pol_coefficients_n *= (pol_degrees[i] + 1); // +1 because starts from 0
	}

	pol_coefficients_n -= 1;                      // because of a(0...0) = 1

	// Initialize conditions and marginals

	if (given_ind.size() == 0)                    // if there is no conditioned components
	{
		given_ind = LogicalVector(pol_degrees_n);   // false by default
	}

	if (omit_ind.size() == 0)                     // if there is no marginalized components
	{
		omit_ind = LogicalVector(pol_degrees_n);    // false by default
	}

	// Determine whether initial values have been manually provided

	bool x0_given = true; 

	if (x0.size() == 0) // if x0 has not been provided manually
	{
		x0_given = false;
	  
		x0 = NumericVector(pol_coefficients_n + // 2 * pol_degrees_n since every random vector components
		                   2 * pol_degrees_n);  // introduces additional mean and sd parameters pair
																	         
	}
	
	// Initialize additional variable which helps
	// to assign indices of parameters in x0
	
	int k = 0;

	// Assign indexes for

		// polynomial coefficients
	
	NumericVector pol_coefficients_ind(pol_coefficients_n);

	for (int i = 0; i < pol_coefficients_n; i++)
	{
	  if (!x0_given) // if user has not provided x0 manually then set mean parameter to sample mean
	  {
	    x0[i] = 0;
	  }
	  
		pol_coefficients_ind[i] = i;
	}

		// mean vector

	NumericVector mean_ind(pol_degrees_n);

	for (int i = pol_coefficients_n; i < (pol_coefficients_n + pol_degrees_n); i++)
	{
		mean_ind[k] = i;
	  
		if (!x0_given) // if user has not provided x0 manually then 
		               // set mean parameter to sample mean
		{
			x0[i] = mean(x(_, k));
		}
		
		k++;
	}
	
		// sd vector

	k = 0;

	NumericVector sd_ind(pol_degrees_n);

	for (int i = (pol_coefficients_n + pol_degrees_n); 
       i < (pol_coefficients_n + 2 * pol_degrees_n); i++)
	{
		sd_ind[k] = i;
	  
		if (!x0_given) // if user has not provided x0 manually then set sd parameter 
		               // to sample standard deviation
		{
			x0[i] = sd(x(_, k));
		}
		
		k++;
	}

	// Deal with truncation

	NumericMatrix tr_left_mat(1, pol_degrees_n);
	NumericMatrix tr_right_mat(1, pol_degrees_n);

	if ((tr_left.size() > 0) | (tr_right.size() > 0))
	{
		if (tr_left.size() == 0) // if there is no left truncation set it to negative infinity
		{
			tr_left = NumericVector(pol_degrees_n, R_NegInf);
		}

		if (tr_right.size() == 0) // if there is no right truncation set it to infinity
		{
			tr_right = NumericVector(pol_degrees_n, R_PosInf);
		}

		tr_left_mat(0, _) = tr_left;
		tr_right_mat(0, _) = tr_right;
	} 
	else
	{
		std::fill(tr_left_mat.begin(), tr_left_mat.end(), NA_REAL);
		std::fill(tr_right_mat.begin(), tr_right_mat.end(), NA_REAL);
	}
	
	// Apply optimization routine
	
	  // Set optim control parameters

	List PGN_control = List::create(
	     Named("maxit") = 100000000, 
       Named("fnscale") = -1.0,
       Named("abstol") = std::sqrt(std::numeric_limits<double>::epsilon()),
       Named("reltol") = std::sqrt(std::numeric_limits<double>::epsilon())
	  );
	
	  // Perform the optimization

	List optim_results = optim(
	    Rcpp::_["par"] = x0,
	    Rcpp::_["fn"] = Rcpp::InternalFunction(&hpaLnLOptim),
	    Rcpp::_["gr"] = Rcpp::InternalFunction(&hpaLnLOptim_grad),
	    Rcpp::_["control"] = PGN_control,
	    Rcpp::_["method"] = "BFGS",
	    Rcpp::_["hessian"] = true,
	    Rcpp::_["x_data"] = x,
	    Rcpp::_["pol_coefficients_ind"] = pol_coefficients_ind,
	    Rcpp::_["pol_degrees"] = pol_degrees,
	    Rcpp::_["given_ind"] = given_ind,
	    Rcpp::_["omit_ind"] = omit_ind,
	    Rcpp::_["mean_ind"] = mean_ind,
	    Rcpp::_["sd_ind"] = sd_ind,
	    Rcpp::_["tr_left"] = tr_left_mat,
	    Rcpp::_["tr_right"] = tr_right_mat,
	    Rcpp::_["is_parallel"] = is_parallel);
	
	// Extract optimization results and assign them to the variables
	// representing estimated parameters

	NumericVector x1 = optim_results["par"];
	
	int x1_n = x1.size();

	double lnL = optim_results["value"];

	NumericVector mean = x1[mean_ind];
	NumericVector sd = x1[sd_ind];

	NumericVector pol_coefficients = x1[pol_coefficients_ind];
	pol_coefficients.push_front(1);

	// Get covariance matrix estimate of "cov_type" type
	
	NumericMatrix cov_mat;  // covariance matrix
	
	arma::mat H_part;       // butter of sandwich estimator
	
	arma::mat J_part;       // bread of sandwich estimator
	
	  // Estimate jacobian for the inner part
	
	if ((cov_type == "gop") | (cov_type == "sandwich") | (cov_type == "sandwichSR"))
	{
	  NumericMatrix my_jacobian = hpaLnLOptim_grad_ind(
                      	                   x1, x,
                                           pol_coefficients_ind,
                                           pol_degrees,
                                           given_ind, omit_ind,
                                           mean_ind, sd_ind,
                                           tr_left_mat, tr_right_mat,
                                           is_parallel);
	  
	  J_part = as<arma::mat>(my_jacobian);
	}

	if ((cov_type == "hessian") | (cov_type == "sandwich"))
	{
	  NumericMatrix my_hessian = optim_results["hessian"];

	  H_part = as<arma::mat>(my_hessian).i();
	}
	
	if (cov_type == "sandwich")
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
	  
	    // list to store bootstrap results
	  
	  List boot_List;
	  
	    // bootstrap procedure
	  
	 if (cov_type == "bootstrap")
	 {
	   for(int i = 0; i < boot_iter; i++)
	   {
	     // Generate sample with replacement
	     
	     NumericVector sample_ind = floor(runif(n_obs, 0, n_obs));
	     
	     NumericMatrix boot_sample = NumericMatrix(n_obs, pol_degrees_n);
	     
	     for (int j = 0; j < n_obs; j++)
	     {
	       boot_sample(j, _) = x(sample_ind[j], _);
	     }
	     
	     // Perform estimaton
	     
	     List boot_results = optim(
	       Rcpp::_["par"] = x1,
	       Rcpp::_["fn"] = Rcpp::InternalFunction(&hpaLnLOptim),
	       Rcpp::_["gr"] = Rcpp::InternalFunction(&hpaLnLOptim_grad),
	       Rcpp::_["control"] = PGN_control,
	       Rcpp::_["method"] = "BFGS",
	       Rcpp::_["hessian"] = false,
	       Rcpp::_["x_data"] = boot_sample,
	       Rcpp::_["pol_coefficients_ind"] = pol_coefficients_ind,
	       Rcpp::_["pol_degrees"] = pol_degrees,
	       Rcpp::_["given_ind"] = given_ind,
	       Rcpp::_["omit_ind"] = omit_ind,
	       Rcpp::_["mean_ind"] = mean_ind,
	       Rcpp::_["sd_ind"] = sd_ind,
	       Rcpp::_["tr_left"] = tr_left_mat,
	       Rcpp::_["tr_right"] = tr_right_mat,
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

	NumericMatrix results(x1_n, 3);

	StringVector results_cols = StringVector::create(
	  "Estimate", "Std. Error", "P(>|z|)");
	StringVector results_rows(x1_n);

		// get vector index matrix for polynomial coefficients

	NumericMatrix pol_ind = polynomialIndex(pol_degrees);

		// assign results matrix columns with values for

			// polynomial coefficients

	for (int i = 1; i < (pol_coefficients_n + 1); i++)
	{
		results_rows[(i - 1)] = "a";
		for (int j = 0; j < pol_degrees_n; j++)
		{
			int my_int = pol_ind(j, i);
			results_rows[(i - 1)] = as<std::string>(results_rows[(i - 1)]) + 
			                                        "_" + std::to_string(my_int);
		}
		results((i - 1), 0) = pol_coefficients[i];
		results((i - 1), 1) = sqrt(cov_mat((i - 1), (i - 1)));
		double t_stat = results((i - 1), 0) / results((i - 1), 1);
		NumericVector F_t_stat = pnorm(NumericVector::create(t_stat));
		results((i - 1), 2) = 2 * std::min(F_t_stat[0], 1 - F_t_stat[0]);
	}

			// mean

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

			//	sd

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

		// assign names to rows and some vectors

	rownames(results) = results_rows;
	colnames(results) = results_cols;

	x1.names() = results_rows;

	rownames(cov_mat) = results_rows;
	colnames(cov_mat) = results_rows;

	pol_coefficients.names() = c_R("a_0", results_rows[pol_coefficients_ind]);
	
	if (cov_type == "bootstrap")
	{
  	colnames(boot_parameters) = results_rows;
  	
  	sd_dev.names() = results_rows;
	}

	// Collect the results

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
		Named("n_obs") = n_obs,
		Named("bootstrap") = boot_List);

	// Assign the class to the output list
	
	return_result.attr("class") = "hpaML";

	return(return_result);
}

// Perform log-likelihood function estimation for 
// Phillips-Gallant-Nychka distribution at point
List hpaLnLOptim_List(NumericVector x0,
	NumericMatrix x_data = NumericMatrix(1, 1),
	NumericVector pol_coefficients_ind = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector mean_ind = NumericVector(0),
	NumericVector sd_ind = NumericVector(0),
	NumericMatrix tr_left = NumericMatrix(1, 1),
	NumericMatrix tr_right = NumericMatrix(1, 1),
	bool is_parallel = false)
{
	// Assign values based on their indecies
	
	NumericVector mean = x0[mean_ind];

	NumericVector sd = x0[sd_ind];

	NumericVector pol_coefficients = x0[pol_coefficients_ind];
	pol_coefficients.push_front(1);

	// Initialize values to return 
	
	List return_List;
	
	NumericVector return_individual;
	
	double return_aggregate;

	// Perform calculations

	if (!(R_IsNA(tr_left(0, 0))) & !(R_IsNA(tr_right(0, 0))))
	{
	  return_individual = log(dtrhpa(x_data,
                                   tr_left, tr_right,
                                   pol_coefficients, pol_degrees,
                                   given_ind, omit_ind,
                                   mean, sd,
                                   is_parallel));
	  
	  return_aggregate = sum(return_individual);
	  
	  // if there is some problems related to precision when truncation
	  // has been incorporated then return great negative number
	  
	  if(ihpa(tr_left, tr_right,
	    pol_coefficients, pol_degrees,
	    given_ind, omit_ind, 
	    mean, sd,
	    is_parallel)[0] < std::sqrt(std::numeric_limits<double>::epsilon()))
	  {
	    std::fill(return_individual.begin(), return_individual.end(), -(1e+100));
	    
	    return_List = List::create(Named("aggregate") = -(1e+100),
                                 Named("individual") = return_individual);
	    
	    return(return_List);
	  }
	  
	  return_List = List::create(Named("aggregate") = return_aggregate,
                               Named("individual") = return_individual);

		return(return_List);
	}
	
	return_individual = log(dhpa(x_data,
                  		         pol_coefficients, pol_degrees,
                  		         given_ind, omit_ind,
                  		         mean, sd, 
                  		         is_parallel));
	
	return_aggregate = sum(return_individual);
	
	return_List = List::create(Named("aggregate") = return_aggregate,
                             Named("individual") = return_individual);

	return(return_List);
}

// Get aggregate component (log-lkelihood function value) from hpaLnLOptim_List
double hpaLnLOptim(NumericVector x0,
                   NumericMatrix x_data = NumericMatrix(1, 1),
                   NumericVector pol_coefficients_ind = NumericVector(0),
                   NumericVector pol_degrees = NumericVector(0),
                   LogicalVector given_ind = LogicalVector(0),
                   LogicalVector omit_ind = LogicalVector(0),
                   NumericVector mean_ind = NumericVector(0),
                   NumericVector sd_ind = NumericVector(0),
                   NumericMatrix tr_left = NumericMatrix(1, 1),
                   NumericMatrix tr_right = NumericMatrix(1, 1),
                   bool is_parallel = false)
{
  List return_List = hpaLnLOptim_List(x0, x_data,
                                      pol_coefficients_ind,
                                      pol_degrees,
                                      given_ind, omit_ind,
                                      mean_ind, sd_ind,
                                      tr_left, tr_right,
                                      is_parallel);
  
  double return_aggregate = return_List["aggregate"];
  
  return(return_aggregate);
}

// Get individual component (log-lkelihood function contributions) from hpaLnLOptim_List
NumericVector hpaLnLOptim_ind(NumericVector x0,
                   NumericMatrix x_data = NumericMatrix(1, 1),
                   NumericVector pol_coefficients_ind = NumericVector(0),
                   NumericVector pol_degrees = NumericVector(0),
                   LogicalVector given_ind = LogicalVector(0),
                   LogicalVector omit_ind = LogicalVector(0),
                   NumericVector mean_ind = NumericVector(0),
                   NumericVector sd_ind = NumericVector(0),
                   NumericMatrix tr_left = NumericMatrix(1, 1),
                   NumericMatrix tr_right = NumericMatrix(1, 1),
                   bool is_parallel = false)
{
  List return_List = hpaLnLOptim_List(x0, x_data,
                                      pol_coefficients_ind,
                                      pol_degrees,
                                      given_ind, omit_ind,
                                      mean_ind, sd_ind,
                                      tr_left, tr_right,
                                      is_parallel);
  
  NumericVector return_individual = return_List["individual"];
  
  return(return_individual);
}

// Perform log-likelihood function gradient estimation 
// for Phillips-Gallant-Nychka distribution at point
List hpaLnLOptim_grad_List(NumericVector x0,
                   NumericMatrix x_data = NumericMatrix(1, 1),
                   NumericVector pol_coefficients_ind = NumericVector(0),
                   NumericVector pol_degrees = NumericVector(0),
                   LogicalVector given_ind = LogicalVector(0),
                   LogicalVector omit_ind = LogicalVector(0),
                   NumericVector mean_ind = NumericVector(0),
                   NumericVector sd_ind = NumericVector(0),
                   NumericMatrix tr_left = NumericMatrix(1, 1),
                   NumericMatrix tr_right = NumericMatrix(1, 1),
                   bool is_parallel = false)
{
  
  List return_List;
  
  // Get parameters number
  
  int n_param = x0.size();
  
  int n_obs = x_data.nrow();
  
  // Initialize vector to store gradient (jacobian) values
  
  NumericMatrix my_grad = NumericMatrix(n_obs, n_param);
  
  // Assign values based on their indecies
  
  NumericVector mean = x0[mean_ind];
  
  NumericVector sd = x0[sd_ind];
  
  NumericVector pol_coefficients = x0[pol_coefficients_ind];
  pol_coefficients.push_front(1);
  
  int pol_coefficients_n = pol_coefficients.size();
  
  // Initial function value 
  
  NumericVector fn_values = dhpa(x_data,
                                 pol_coefficients, pol_degrees,
                                 given_ind, omit_ind,
                                 mean, sd,
                                 is_parallel);
  
  // Denominator for truncated densities
  
  NumericVector fn_values_cdf;
  
  if (!(R_IsNA(tr_left(0, 0))) & !(R_IsNA(tr_right(0, 0))))
  {
    fn_values_cdf = ihpa(tr_left, tr_right,
                         pol_coefficients, pol_degrees,
                         given_ind, omit_ind,
                         mean, sd,
                         is_parallel);
    
    if (fn_values_cdf[0] < std::sqrt(std::numeric_limits<double>::epsilon()))
    {
      std::fill(my_grad.begin(), my_grad.end(),
                -(1e+100));
      
      return_List = List::create(Named("aggregate") = colSums(my_grad),
                                 Named("individual") = my_grad);
      
      return(return_List);
    }
  }
  
  // Analytical part of gradient (for polynomial coefficients)
  
  NumericMatrix pol_coefficients_grad = dhpaDiff(x_data,
                                                 pol_coefficients, pol_degrees,
                                                 given_ind, omit_ind,
                                                 mean, sd,
                                                 "pol_coefficients",
                                                 is_parallel);
  
  for (int i = 0; i < (pol_coefficients_n - 1); i++) // for each parameter
  {
    my_grad(_, i) = pol_coefficients_grad(_, i + 1) / fn_values;
  }
  
    // Deal with truncation
  
  NumericMatrix pol_coefficients_grad_cdf;
  
  if (!(R_IsNA(tr_left(0, 0))) & !(R_IsNA(tr_right(0, 0))))
  {
    pol_coefficients_grad_cdf = ihpaDiff(tr_left, tr_right,
                                         pol_coefficients, pol_degrees,
                                         given_ind, omit_ind,
                                         mean, sd,
                                         "pol_coefficients",
                                         is_parallel);
    
    for (int i = 0; i < (pol_coefficients_n - 1); i++) // for each parameter
    {
      double tr_grad = sum(pol_coefficients_grad_cdf(_, i + 1) / fn_values_cdf);
      my_grad(_, i) = my_grad(_, i) - tr_grad;
    }
  }
  
  // Numeric part of gradient (for mu and sigma)
  
    // Set differentiation increment value for each parameter
  
  double machinePrecision = std::numeric_limits<double>::epsilon();
  double my_precision = std::sqrt(machinePrecision);
  
  NumericVector eps = abs(x0 * my_precision);
  
    // Control for zero values
  
  eps[eps < (machinePrecision * 100)] = my_precision;
  
    // Estimate the gradient itself
  
  NumericVector x0_eps = clone(x0);
  
  NumericVector f1;
  
  NumericVector f2;
  
  for (int i = pol_coefficients_n - 1; i < n_param; i++) // for each parameter
  {
      
      // Calculate f(x - eps)
      
      x0_eps[i] = x0[i] - eps[i];

      f1 = hpaLnLOptim_ind(x0_eps, x_data,
                           pol_coefficients_ind,
                           pol_degrees,
                           given_ind, omit_ind,
                           mean_ind, sd_ind = sd_ind,
                           tr_left, tr_right,
                           is_parallel
        );
    
      // Calculate f(x + eps)
    
      x0_eps[i] = x0[i] + eps[i];

      f2 = hpaLnLOptim_ind(x0_eps, x_data,
                           pol_coefficients_ind,
                           pol_degrees,
                           given_ind, omit_ind,
                           mean_ind, sd_ind,
                           tr_left, tr_right,
                           is_parallel
        );
      
      // Estimate the gradient 
        
      my_grad(_, i) = (f2 - f1) / (2 * eps[i]);
    
      // Set x0_eps value to default
    
      x0_eps[i] = x0[i];
  }
  
  return_List = List::create(Named("aggregate") = colSums(my_grad),
                             Named("individual") = my_grad);

  return(return_List);
}

// Get aggregaate component (log-lkelihood function gradient) 
// from hpaLnLOptim_grad_List
NumericVector hpaLnLOptim_grad(NumericVector x0,
                           NumericMatrix x_data = NumericMatrix(1, 1),
                           NumericVector pol_coefficients_ind = NumericVector(0),
                           NumericVector pol_degrees = NumericVector(0),
                           LogicalVector given_ind = LogicalVector(0),
                           LogicalVector omit_ind = LogicalVector(0),
                           NumericVector mean_ind = NumericVector(0),
                           NumericVector sd_ind = NumericVector(0),
                           NumericMatrix tr_left = NumericMatrix(1, 1),
                           NumericMatrix tr_right = NumericMatrix(1, 1),
                           bool is_parallel = false)
{
  List return_List = hpaLnLOptim_grad_List(x0, x_data,
                                      pol_coefficients_ind,
                                      pol_degrees,
                                      given_ind, omit_ind,
                                      mean_ind, sd_ind,
                                      tr_left, tr_right,
                                      is_parallel);
  
  NumericVector return_aggregate = return_List["aggregate"];
  
  return(return_aggregate);
}

// Get individual component (log-lkelihood function gradient contributions) 
// from hpaLnLOptim_grad_List
NumericMatrix hpaLnLOptim_grad_ind(NumericVector x0,
                               NumericMatrix x_data = NumericMatrix(1, 1),
                               NumericVector pol_coefficients_ind = NumericVector(0),
                               NumericVector pol_degrees = NumericVector(0),
                               LogicalVector given_ind = LogicalVector(0),
                               LogicalVector omit_ind = LogicalVector(0),
                               NumericVector mean_ind = NumericVector(0),
                               NumericVector sd_ind = NumericVector(0),
                               NumericMatrix tr_left = NumericMatrix(1, 1),
                               NumericMatrix tr_right = NumericMatrix(1, 1),
                               bool is_parallel = false)
{
  List return_List = hpaLnLOptim_grad_List(x0, x_data,
                                           pol_coefficients_ind,
                                           pol_degrees,
                                           given_ind, omit_ind,
                                           mean_ind, sd_ind,
                                           tr_left, tr_right,
                                           is_parallel);
  
  NumericMatrix return_individual = return_List["individual"];
  
  return(return_individual);
}

//' Predict method for hpaML
//' @param object Object of class "hpaML"
//' @template newdata_Template
//' @return This function returns predictions based 
//' on \code{\link[hpa]{hpaML}} estimation results.
//' @export
// [[Rcpp::export]]
NumericVector predict_hpaML(List object, 
                            NumericMatrix newdata = NumericMatrix(1, 1))
{
	List model = object;

	// Load additional environments

		// stats environment
		
	Rcpp::Environment stats_env("package:stats");
	Rcpp::Function na_omit_R = stats_env["na.omit"];

	// Get distribution parameters

	NumericVector pol_coefficients = model["pol_coefficients"];
	NumericVector pol_degrees = model["pol_degrees"];

	NumericVector mean = model["mean"];
	NumericVector sd = model["sd"];

	NumericMatrix tr_left = model["tr_left"];
	NumericMatrix tr_right = model["tr_right"];

	LogicalVector omit_ind = model["omit_ind"];
	LogicalVector given_ind = model["given_ind"];

	NumericMatrix x = model["data"];

	// Get data
	
	if ((newdata.ncol() == 1) & (newdata.nrow() == 1))
	{
		newdata = x;
	} 
	else 
	{
		newdata = na_omit_R(newdata);
	}

	// Estimate
	
	if (!(R_IsNA(tr_left(0, 0))) & !(R_IsNA(tr_right(0, 0))))
	{
		return(dtrhpa(newdata,
			tr_left, tr_right,
			pol_coefficients, pol_degrees,
			given_ind, omit_ind,
			mean, sd, false));
	}

	return(dhpa(newdata,
		pol_coefficients, pol_degrees,
		given_ind, omit_ind,
		mean, sd, false));
}

//' Summarizing hpaML Fits
//' @param object Object of class "hpaML"
//' @return This function returns the same 
//' list as \code{\link[hpa]{hpaML}} function changing 
//' it's class to "summary.hpaML".
//' @export
// [[Rcpp::export]]
List summary_hpaML(List object)
{
	List return_result = clone(object); // in order to preserve model class

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

	// Load additional environments

	Rcpp::Environment base_env("package:base");
	Rcpp::Function as_table = base_env["as.table"];
	Rcpp::Function cbind = base_env["cbind"];
	Rcpp::Function round_R = base_env["round"];
	Rcpp::Function print_R = base_env["print"];
	Rcpp::Function cat_R = base_env["cat"];

	// Do other obvious stuff

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
	std::string df_string = std::to_string(df) + 
	  " free parameters (df = " + std::to_string(n_obs - df) + ")" + "\n";

	cat_R("--------------------------------------------------------------\n");
	
	cat_R("Semi-nonparametric maximum likelihood estimation\n");
	
	cat_R("---\n");

	cat_R(lnL_string);
	cat_R(AIC_string);
	cat_R(n_obs_string);
	cat_R(df_string);
	
	cat_R("---\n");

	cat_R("Distribution parameters:\n");
	print_R(as_table(cbind(results, stars)));

	cat_R("---\n");
	cat_R("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n");

	cat_R("--------------------------------------------------------------\n");
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

//' Calculates multivariate empirical cumulative distribution function
//' @description This function calculates multivariate 
//' empirical cumulative distribution function
//' at each point of the sample
//' @param x numeric matrix which rows are observations
//' @export
// [[Rcpp::export]]
NumericVector mecdf(NumericMatrix x)
{
  int n = x.nrow();
  int m = x.ncol();
  
  NumericVector my_cdf = NumericVector(n);
  
  for (int i = 0; i < n; i++)
  {
    for (int j = i; j < n; j++)
    {
      
      int n_greater = 0;
      
      for (int t = 0; t < m; t++)
      {
        if (x(j, t) <= x(i, t))
        {
          n_greater++;
        }
      }
      
      if (n_greater == m)
      {
        my_cdf[i]++;
      }
      
      if (n_greater == 0)
      {
        my_cdf[j]++;
      }
      
    }
  }
  
  my_cdf = my_cdf / n;
  
  return(my_cdf);
}
