#include "normalMoments.h"
#include "polynomialIndex.h"
#include "hpaMain.h"
#include "ParallelFunctions.h"
#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace RcppParallel;

// Hermite polynomial density,
// cumulative distribution function and moments approximations.
// @description This function estimates hermite polynomial density,
// cumulative distribution function and moments approximations.
// @template x_lower_Template
// @template x_upper_Template
// @template pol_coefficients_Template
// @template pol_degrees_Template
// @template type_Template
// @template given_ind_Template
// @template omit_ind_Template
// @template mean_Template
// @template sd_Template
// @template expectation_powers_Template
// @details
// If you already have some precalculated values please specify them using 
// \code{pdf_lower}, \code{pdf_upper}, \code{cdf_lower}, \code{cdf_upper} and \code{cdf_difference} arguments.
List hpaMain(
	NumericMatrix x_lower = NumericMatrix(1,1),
	NumericMatrix x_upper = NumericMatrix(1,1),
	NumericVector pol_coefficients = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	String type = "pdf",
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector mean = NumericVector(0),
	NumericVector sd = NumericVector(0),
	NumericVector expectation_powers = NumericVector(0),
	String grad_type = "NO",
	bool is_parallel = false,
	bool is_cdf = false)
{

	// Get number of observations
	int n = x_upper.nrow();

	// Initialize polynomial structure related values
	int pol_degrees_n = pol_degrees.size();
	int pol_coefficients_n = pol_coefficients.size();

	// Fill x_lower with (-INF) if need
	if (!((type == "interval") | (type == "expectation truncated")) | is_cdf)
	{
		x_lower = NumericMatrix(n, pol_degrees_n);
		std::fill(x_lower.begin(), x_lower.end(), R_NegInf);
	}

	// Initialize conditions and marginals
	if (given_ind.size() == 0)
	{
		given_ind = LogicalVector(pol_degrees_n);
	}
	if (omit_ind.size() == 0)
	{
		omit_ind = LogicalVector(pol_degrees_n);
	}
	
	// Initialize mean and standard deviations
	if (mean.size() == 0)
	{
		mean = NumericVector(pol_degrees_n);
		std::fill(mean.begin(), mean.end(), 0);
	}

	if (sd.size() == 0)
	{
		sd = NumericVector(pol_degrees_n);
		std::fill(sd.begin(), sd.end(), 1);
	}

	// Control for the expected powered product powers values
	if ((expectation_powers.size() == 0) | ((type != "expectation") & (type != "expectation truncated")))
	{
		expectation_powers = NumericVector(pol_degrees_n);
		std::fill(expectation_powers.begin(), expectation_powers.end(), 0);
	}

	// Initialize indexes for observable unconditioned components
	LogicalVector d_cond = ((!given_ind) & (!omit_ind));
	
	// Define vectors related to normal distribution
	
	  // pdf associated values
	NumericMatrix pdf_upper = NumericMatrix(n, pol_degrees_n);
	NumericMatrix pdf_lower = NumericMatrix(n, pol_degrees_n);
	NumericVector pdf_product(n);
	  
	  // cdf associated values
	NumericMatrix cdf_upper = NumericMatrix(n, pol_degrees_n);
	NumericMatrix cdf_lower = NumericMatrix(n, pol_degrees_n);
	NumericMatrix cdf_difference = NumericMatrix(n, pol_degrees_n);
	NumericVector cdf_difference_product = NumericVector(n);
	std::fill(cdf_difference_product.begin(), cdf_difference_product.end(), 1);

	if (type != "expectation")
	{
		// Initialize densities

  		// Upper densities
  		if (pdf_upper(0, 0) == 0)
  		{
  			pdf_upper = NumericMatrix(n, pol_degrees_n);
  			for (int i = 0; i < pol_degrees_n; i++)
  			{
  				if (!omit_ind[i])
  				{
  				  pdf_upper(_, i) = dnorm_parallel(x_upper(_, i), mean[i], sd[i], is_parallel);
  				}
  			}
  		}
  
  		// Lower densities
  		if ((type == "interval") | (type == "expectation truncated"))
  		{
  			for (int i = 0; i < pol_degrees_n; i++)
  			{
  				if (!omit_ind[i])
  				{
  				  pdf_lower(_, i) = dnorm_parallel(x_lower(_, i), mean[i], sd[i], is_parallel);
  				}
  			}
  		}

  		// Product of densities
  		if (type == "pdf")
  		{
  			std::fill(pdf_product.begin(), pdf_product.end(), 1);
  			for (int i = 0; i < pol_degrees_n; i++)
  			{
  				if (d_cond[i])
  				{
  					pdf_product = pdf_product * pdf_upper(_, i);
  				}
  			}
  		}

		// Initialize cumulative distribution functions (cdfs)
		
  		// Upper cdf
  		if (type != "pdf")
  		{
  			if (cdf_upper(0, 0) == 0)
  			{
  				for (int i = 0; i < pol_degrees_n; i++)
  				{
  					if (d_cond[i])
  					{
  					  cdf_upper(_, i) = pnorm_parallel(x_upper(_, i), mean[i], sd[i], is_parallel);
  					}
  				}
  			}
  		}

  			// Lower cdf
  			if (((type == "interval") | (type == "expectation truncated")))
  			{
  				cdf_lower = NumericMatrix(n, pol_degrees_n);
  				for (int i = 0; i < pol_degrees_n; i++)
  				{
  					if (d_cond[i])
  					{
  					  cdf_lower(_, i) = pnorm_parallel(x_lower(_, i), mean[i], sd[i], is_parallel);
  					}
  				}
  			}

  			// Calculate cdf_difference
  			for (int i = 0; i < pol_degrees_n; i++)
  			{
  				if (d_cond[i])
  				{
  					cdf_difference(_, i) = (cdf_upper(_, i) - cdf_lower(_, i));
  				}
  			}

  		// Estimate cdf_difference product
  		if (type != "pdf")
  		{
  			for (int i = 0; i < pol_degrees_n; i++)
  			{
  				if (d_cond[i])
  				{
  					cdf_difference_product = cdf_difference_product * cdf_difference(_, i);
  				}
  			}
  		}
	}

	// Define vector indexing system for polynomial
	NumericMatrix polynomial_index = polynomialIndex(pol_degrees);

	// Calculate moments
	List moments(pol_degrees_n);
	int max_degree;

	for (int i = 0; i < pol_degrees_n; i++)
	{
		if (!given_ind[i])
		{
			max_degree = 2 * pol_degrees[i] + expectation_powers[i];
			moments[i] = normalMoment(max_degree,
									  mean[i], sd[i], 
									  true, false);
		}
	}

	// Calcule truncated moments
	List tr_moments(pol_degrees_n);

	if ((type != "pdf") & (type != "expectation"))
	{
		for (int i = 0; i < pol_degrees_n; i++)
		{
			if (d_cond[i])
			{
				max_degree = 2 * pol_degrees[i] + expectation_powers[i];
				// Note that arguments order matters!
				tr_moments[i] = truncatedNormalMoment(max_degree,
					x_lower(_, i), x_upper(_, i),
					mean[i], sd[i],
					pdf_lower(_, i), cdf_lower(_, i),
					pdf_upper(_, i), cdf_upper(_, i),
					cdf_difference(_, i), true, false, is_parallel);
			}
		}
	}

	// Calculate x powers
	LogicalVector x_cond(pol_degrees_n);

	if (type == "pdf")
	{
		x_cond = !omit_ind;
	} else {
		x_cond = given_ind;
	}

	List x_pow(pol_degrees_n);

	int k = 0;

	for (int i = 0; i < pol_degrees_n; i++)
	{
		if (x_cond[i])
		{
			k = 2 * pol_degrees[i] + 1;
			x_pow[i] = NumericMatrix(n, k);
			NumericMatrix x_pow_i = x_pow[i]; // it is reference
			NumericVector x_upper_i = x_upper(_, i);
			for (int j = 0; j < k; j++)
			{
			  if(is_parallel)
			  {
			    x_pow_i(_, j) = ParallelVectorPow(x_upper_i, j);
			  } else {
			    x_pow_i(_, j) = pow(x_upper_i, j);
			  }
			}
		}
	}

	// Calculate main expression

	// Initialize values to store temporal results

	NumericVector value_pgn(n);
	std::fill(value_pgn.begin(), value_pgn.end(), 0);

	NumericVector psi(n);
	std::fill(psi.begin(), psi.end(), 0);

	NumericVector value_sum_element(n);
	NumericVector psi_sum_element(n);
	double polynomial_sum = 0;
	
	// Initialize values to store gradient information if need
	
	NumericMatrix pc_grad = NumericMatrix(n, pol_coefficients_n);
	NumericMatrix pc_grad_value = NumericMatrix(n, pol_coefficients_n);
	NumericMatrix pc_grad_psi = NumericMatrix(n, pol_coefficients_n);
	
	if(grad_type == "pol_coefficients")
	{
	  pc_grad = NumericMatrix(n, pol_coefficients_n);
	  pc_grad_value = NumericMatrix(n, pol_coefficients_n);
	  pc_grad_psi = NumericMatrix(n, pol_coefficients_n);
	}

	// Perform main calculations

	for (int i = 0; i < pol_coefficients_n; i++)
	{
		for (int j = i; j < pol_coefficients_n; j++)
		{
		  
		  double pol_coefficients_prod = pol_coefficients[i] * pol_coefficients[j];
		  
			// Initialize temporal value
			std::fill(value_sum_element.begin(), value_sum_element.end(), 1);
			std::fill(psi_sum_element.begin(), psi_sum_element.end(), 1);

			// Main calculations
			for (int r = 0; r < pol_degrees_n; r++)
			{
				polynomial_sum = polynomial_index(r, i) + polynomial_index(r, j);
				if (!omit_ind[r])
				{
					if ((type == "pdf") | (given_ind[r]))
					{
						NumericMatrix x_pow_r = x_pow[r];
						value_sum_element = value_sum_element * x_pow_r(_, polynomial_sum);
					} else {
						if (type != "expectation")
						{
							NumericMatrix tr_moments_r = tr_moments[r];
							value_sum_element = value_sum_element *
								tr_moments_r(_, polynomial_sum + expectation_powers[r]);
						} else {
							NumericVector moments_r_e = moments[r];
							value_sum_element = value_sum_element *
								moments_r_e[polynomial_sum + expectation_powers[r]];
						}
					}
				} else {
					NumericVector moments_r = moments[r];
					value_sum_element = value_sum_element * moments_r[polynomial_sum];
				}
				// psi
				if (given_ind[r])
				{
					NumericMatrix x_pow_r = x_pow[r];
					psi_sum_element = psi_sum_element * x_pow_r(_, polynomial_sum);
				} else {
					if (type != "expectation truncated")
					{
						NumericVector moments_r = moments[r];
						psi_sum_element = psi_sum_element * moments_r[polynomial_sum];
					} else {
						NumericMatrix tr_moments_r = tr_moments[r];
						psi_sum_element = psi_sum_element * tr_moments_r(_, polynomial_sum);
					}
				}
			}
			// Each iteration perform results storage
			int mult_for_unequal_i_j = (1 + (i != j));
			value_pgn = value_pgn + mult_for_unequal_i_j * value_sum_element * pol_coefficients_prod;
			psi = psi + mult_for_unequal_i_j * psi_sum_element * pol_coefficients_prod;
			
			// gradient specific storage if need
			if(grad_type == "pol_coefficients")
			{
			  pc_grad_value(_, i) = pc_grad_value(_, i) + mult_for_unequal_i_j * value_sum_element * pol_coefficients[j];
			  pc_grad_value(_, j) = pc_grad_value(_, j) + mult_for_unequal_i_j * value_sum_element * pol_coefficients[i];
			  pc_grad_psi(_, i) = pc_grad_psi(_, i) + mult_for_unequal_i_j * psi_sum_element * pol_coefficients[j];
			  pc_grad_psi(_, j) = pc_grad_psi(_, j) + mult_for_unequal_i_j * psi_sum_element * pol_coefficients[i];
			}
		}
	}

	// Return gradient specific values if need

	if (grad_type == "pol_coefficients")
	{
	  NumericVector psi_2 = psi * psi;
	  
	    if(type == "pdf")
	    {
	      NumericVector value_2 = value_pgn * pdf_product;
	      for (int i = 0; i < pol_coefficients_n; i++)
	      {
	        pc_grad(_, i) = ((pc_grad_value(_, i) * pdf_product) / psi) -
	          (pc_grad_psi(_, i) * value_2 / psi_2);
	      }
	      return(List::create(Named("grads") = pc_grad));
	    }
	    
	    if(type == "interval")
	    {
	      NumericVector value_2 = value_pgn * cdf_difference_product;
	      for (int i = 0; i < pol_coefficients_n; i++)
	      {
	        pc_grad(_, i) = ((pc_grad_value(_, i) * cdf_difference_product) / psi) -
	          (pc_grad_psi(_, i) * value_2 / psi_2);
	      }
	      return(List::create(Named("grads") = pc_grad));
	    }
	}

	// Return the result depending on the type of calculations

	if ((type == "expectation") | (type == "expectation truncated"))
	{
		return(List::create(Named("values") = value_pgn / psi));
	}

	if (type != "pdf")
	{
	  return(List::create(Named("values") = (value_pgn * cdf_difference_product) / psi));
	} else {
	  return(List::create(Named("values") = (value_pgn * pdf_product) / psi));
	}

	return(List::create(Named("values") = NumericVector(0)));
}

//' Density function hermite polynomial approximation
//' @description This function calculates density function hermite polynomial approximation.
//' @template x_pdf_Template
//' @template pol_coefficients_Template
//' @template pol_degrees_Template
//' @template given_ind_Template
//' @template omit_ind_Template
//' @template mean_Template
//' @template sd_Template
//' @template is_parallel_Template
//' @template GN_details_Template
//' @template dhpa_examples_Template
//' @return This function returns density function hermite polynomial approximation at point \code{x}.
//' @export
// [[Rcpp::export]]
NumericVector dhpa(
	NumericMatrix x = NumericMatrix(1,1),
	NumericVector pol_coefficients = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector mean = NumericVector(0),
	NumericVector sd = NumericVector(0),
	bool is_parallel = false)
{

  List return_List = hpaMain(
    NumericMatrix(1, 1),            // x_lower
    x,                              // x_upper
    pol_coefficients, pol_degrees,
    "pdf",                          // type
    given_ind, omit_ind,
    mean, sd,
    NumericVector(0),"NO", 
    is_parallel);
		                   
		NumericVector return_value = return_List["values"];
		                   
		return(return_value);
}

//' Distribution function hermite polynomial approximation
//' @description This function calculates cumulative distribution function hermite polynomial approximation.
//' @template x_cdf_Template
//' @template pol_coefficients_Template
//' @template pol_degrees_Template
//' @template given_ind_Template
//' @template omit_ind_Template
//' @template mean_Template
//' @template sd_Template
//' @template is_parallel_Template
//' @template GN_details_Template
//' @return This function returns cumulative distribution function hermite polynomial approximation at point \code{x}.
//' @template phpa_examples_Template
//' @export
// [[Rcpp::export]]
NumericVector phpa(
	NumericMatrix x = NumericMatrix(1,1),
	NumericVector pol_coefficients = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector mean = NumericVector(0),
	NumericVector sd = NumericVector(0),
	bool is_parallel = false) 
{
  
  List return_List = hpaMain(
    NumericMatrix(1, 1),                     // x_lower
    x,                                       // x_upper
    pol_coefficients, pol_degrees,
    "cdf",                                   // type
    given_ind, omit_ind,
    mean, sd,
    NumericVector(0),"NO", 
    is_parallel, true);
  
  NumericVector return_value = return_List["values"];
  
  return(return_value);
}

//' Interval distribution function hermite polynomial approximation
//' @description This function calculates interval distribution function hermite polynomial approximation.
//' @template x_lower_Template
//' @template x_upper_Template
//' @template pol_coefficients_Template
//' @template pol_degrees_Template
//' @template given_ind_Template
//' @template omit_ind_Template
//' @template mean_Template
//' @template sd_Template
//' @template is_parallel_Template
//' @template interval_cdf_Template
//' @template GN_details_Template
//' @return This function returns interval distribution function hermite polynomial approximation at point \code{x}.
//' @template ihpa_examples_Template
//' @export
// [[Rcpp::export]]
NumericVector ihpa(
	NumericMatrix x_lower = NumericMatrix(1, 1),
	NumericMatrix x_upper = NumericMatrix(1, 1),
	NumericVector pol_coefficients = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector mean = NumericVector(0),
	NumericVector sd = NumericVector(0),
	bool is_parallel = false) 
{

  List return_List = hpaMain(
    x_lower,                                 // x_lower
    x_upper,                                 // x_upper
    pol_coefficients, pol_degrees,
    "interval",                              // type
    given_ind, omit_ind,
    mean, sd,
    NumericVector(0),"NO", 
    is_parallel);
  
  NumericVector return_value = return_List["values"];
  
  return(return_value);
}

//' Expected powered product hermite polynomial approximation
//' @description This function calculates expected powered product hermite polynomial approximation.
//' @template x_expectation_Template
//' @template pol_coefficients_Template
//' @template pol_degrees_Template
//' @template given_ind_Template
//' @template omit_ind_Template
//' @template mean_Template
//' @template sd_Template
//' @template expectation_powers_Template
//' @template is_parallel_Template
//' @template expected_powered_product_Template
//' @template GN_details_Template
//' @return This function returns numeric vector of expected powered product hermite polynomial approximations.
//' @template ehpa_examples_Template
//' @export
// [[Rcpp::export]]
NumericVector ehpa(NumericMatrix x = NumericMatrix(1, 1), //for given
	NumericVector pol_coefficients = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector mean = NumericVector(0),
	NumericVector sd = NumericVector(0),
	NumericVector expectation_powers = NumericVector(0),
	bool is_parallel = false) 
{
  
  List return_List = hpaMain(
    NumericMatrix(1, 1),                     // x_lower
    x,                                       // x_upper
    pol_coefficients, pol_degrees,
    "expectation",                           // type
    given_ind, omit_ind,
    mean, sd,
    expectation_powers,
    "NO", 
    is_parallel);
  
  NumericVector return_value = return_List["values"];
  
  return(return_value);
}

//' Expected powered product hermite polynomial approximation for truncated distribution
//' @description This function calculates expected powered product hermite polynomial approximation for truncated distribution.
//' @template tr_left_Template
//' @template tr_right_Template
//' @template pol_coefficients_Template
//' @template pol_degrees_Template
//' @template mean_Template
//' @template sd_Template
//' @template expectation_powers_Template
//' @template is_parallel_Template
//' @template expected_powered_product_Template
//' @template GN_details_Template
//' @template etrhpa_examples_Template
//' @return This function returns numeric vector of expected powered product hermite polynomial approximations for truncated distribution.
//' @export
// [[Rcpp::export]]
NumericVector etrhpa(
	NumericMatrix tr_left = NumericMatrix(1, 1),
	NumericMatrix tr_right = NumericMatrix(1, 1),
	NumericVector pol_coefficients = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	NumericVector mean = NumericVector(0),
	NumericVector sd = NumericVector(0),
	NumericVector expectation_powers = NumericVector(0),
	bool is_parallel = false) 
{

  List return_List = hpaMain(
    tr_left,                                 // x_lower
    tr_right,                                // x_upper
    pol_coefficients, pol_degrees,
    "expectation truncated",                 // type
    LogicalVector(0), LogicalVector(0),      // given_ind, omit_ind
    mean, sd,
    expectation_powers,
    "NO", 
    is_parallel);
  
  NumericVector return_value = return_List["values"];
  
  return(return_value);
}

//' Truncated density function hermite polynomial approximation
//' @description This function calculates truncated density function hermite polynomial approximation.
//' @template x_pdf_Template
//' @template tr_left_Template
//' @template tr_right_Template
//' @template pol_coefficients_Template
//' @template pol_degrees_Template
//' @template given_ind_Template
//' @template omit_ind_Template
//' @template mean_Template
//' @template sd_Template
//' @template is_parallel_Template
//' @template GN_details_Template
//' @template dtrhpa_examples_Template
//' @return This function returns density function hermite polynomial approximation at point \code{x} for truncated distribution.
//' @export
// [[Rcpp::export]]
NumericVector dtrhpa(
	NumericMatrix x = NumericMatrix(1, 1),
	NumericMatrix tr_left = NumericMatrix(),
	NumericMatrix tr_right = NumericMatrix(),
	NumericVector pol_coefficients = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector mean = NumericVector(0),
	NumericVector sd = NumericVector(0),
	bool is_parallel = false)
{

	// Calculate the nominator
	NumericVector density_main = dhpa(
		x,
		pol_coefficients, pol_degrees,
		given_ind, omit_ind,
		mean, sd,
		is_parallel);

	// Calculate the denonominator
	NumericVector density_tr = ihpa( 
		tr_left, tr_right,
		pol_coefficients, pol_degrees,
		given_ind, omit_ind, 
		mean, sd, 
	  is_parallel);

	if ((tr_left.nrow() == 1) | (tr_right.nrow() == 1))
	{
		return(density_main / density_tr[0]);
	} 

	return(density_main / density_tr);
}

//' Truncated interval distribution function hermite polynomial approximation for truncated distribution
//' @description This function calculates truncated interval distribution function hermite polynomial approximation for truncated distribution.
//' @template x_lower_Template
//' @template x_upper_Template
//' @template tr_left_Template
//' @template tr_right_Template
//' @template pol_coefficients_Template
//' @template pol_degrees_Template
//' @template given_ind_Template
//' @template omit_ind_Template
//' @template mean_Template
//' @template sd_Template
//' @template is_parallel_Template
//' @template itrhpa_examples_Template
//' @template interval_cdf_Template
//' @template GN_details_Template
//' @return This function returns interval distribution function (idf) hermite polynomial approximation at point \code{x} for truncated distribution.
//' @export
// [[Rcpp::export]]
NumericVector itrhpa(
	NumericMatrix x_lower = NumericMatrix(1, 1),
	NumericMatrix x_upper = NumericMatrix(1, 1),
	NumericMatrix tr_left = NumericMatrix(1, 1),
	NumericMatrix tr_right = NumericMatrix(1, 1),
	NumericVector pol_coefficients = NumericVector(0),
	NumericVector pol_degrees = NumericVector(0),
	LogicalVector given_ind = LogicalVector(0),
	LogicalVector omit_ind = LogicalVector(0),
	NumericVector mean = NumericVector(0),
	NumericVector sd = NumericVector(0),
	bool is_parallel = false)
{

	// Calculate the nominator
	NumericVector interval_main = ihpa(
		x_lower, x_upper,
		pol_coefficients, pol_degrees,
		given_ind, omit_ind,
		mean, sd,
    is_parallel);

	// Calculate the denominator
	NumericVector interval_tr = ihpa(
		tr_left, tr_right,
		pol_coefficients, pol_degrees,
		given_ind, omit_ind,
		mean, sd,
		is_parallel);

	if ((tr_left.nrow() == 1) | (tr_right.nrow()))
	{
		return(interval_main / interval_tr[0]);
	}

	return(interval_main / interval_tr);
}

//' Calculate gradient of density function hermite polynomial approximation
//' @description This function calculates gradient of density function 
//' hermite polynomial approximation.
//' @template x_pdf_Template
//' @template pol_coefficients_Template
//' @template pol_degrees_Template
//' @template given_ind_Template
//' @template omit_ind_Template
//' @template mean_Template
//' @template sd_Template
//' @template type_diff_Template
//' @template is_parallel_Template
//' @template GN_details_Template
//' @template dhpaDiff_examples_Template
//' @return This function returns gradient of density function hermite polynomial 
//' approximation at point \code{x}. Gradient elements are determined 
//' by the \code{type} argument.
//' @details
//' If \code{x} has more then one row then the output will be jacobian matrix where 
//' rows are gradients.
//' @export
// [[Rcpp::export]]
NumericMatrix dhpaDiff(
    NumericMatrix x = NumericMatrix(1, 1),
    NumericVector pol_coefficients = NumericVector(0),
    NumericVector pol_degrees = NumericVector(0),
    LogicalVector given_ind = LogicalVector(0),
    LogicalVector omit_ind = LogicalVector(0),
    NumericVector mean = NumericVector(0),
    NumericVector sd = NumericVector(0),
    String type = "pol_coefficients",
    bool is_parallel = false)
{
  
  List return_List;
  
  if (type == "pol_coefficients")
  {
  return_List = hpaMain(
    NumericMatrix(1, 1),            // x_lower
    x,                              // x_upper
    pol_coefficients, pol_degrees,
    "pdf",                          // type
    given_ind, omit_ind,
    mean, sd,
    NumericVector(0),
    type,                           // grad type
    is_parallel);                          
  }
  
  NumericMatrix return_value = return_List["grads"];
  
  return(return_value);
}

//' Calculate gradient of interval distribution function hermite polynomial approximation
//' @description This function calculates interval distribution function hermite polynomial approximation.
//' @template x_lower_Template
//' @template x_upper_Template
//' @template pol_coefficients_Template
//' @template pol_degrees_Template
//' @template given_ind_Template
//' @template omit_ind_Template
//' @template mean_Template
//' @template sd_Template
//' @template type_diff_Template
//' @template is_parallel_Template
//' @template interval_cdf_Template
//' @template GN_details_Template
//' @return This function returns gradient of interval distribution function hermite polynomial 
//' approximation at point \code{x}. Gradient elements are determined 
//' by the \code{type} argument.
//' @details
//' If \code{x} has more then one row then the output will be jacobian matrix where 
//' rows are gradients.
//' @template ihpaDiff_examples_Template
//' @export
// [[Rcpp::export]]
NumericMatrix ihpaDiff(
    NumericMatrix x_lower = NumericMatrix(1, 1),
    NumericMatrix x_upper = NumericMatrix(1, 1),
    NumericVector pol_coefficients = NumericVector(0),
    NumericVector pol_degrees = NumericVector(0),
    LogicalVector given_ind = LogicalVector(0),
    LogicalVector omit_ind = LogicalVector(0),
    NumericVector mean = NumericVector(0),
    NumericVector sd = NumericVector(0),
    String type = "pol_coefficients",
    bool is_parallel = false)
{
  
  List return_List;
  
  if (type == "pol_coefficients")
  {
    return_List = hpaMain(
      x_lower,                        // x_lower
      x_upper,                        // x_upper
      pol_coefficients, pol_degrees,
      "interval",                     // type
      given_ind, omit_ind,
      mean, sd,
      NumericVector(0),
      type,                           // grad type
      is_parallel);                          
  }
  
  NumericMatrix return_value = return_List["grads"];
  
  return(return_value);
}
