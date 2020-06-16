#include "polynomialIndex.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends(RcppArmadillo)]]

//' Returns matrix of polynomial indexes
//' @description Returns matrix of polynomial indexes for the polynomial with degrees (orders) vector \code{pol_degrees}.
//' @template pol_degrees_Template
//' @details
//' This function motivation is to have an opportunity to
//' iterate through the columns of polynomial indexes matrix in order to access polynomial elements 
//' being aware of their powers.
//' @template polynomialIndex_examples_Template
//' @return This function returns polynomial indexes matrix which rows are responsible 
//' for variables while columns are related to powers.
//' @export
// [[Rcpp::export]]
NumericMatrix polynomialIndex(NumericVector pol_degrees = 0) 
{
  // Convert pol_degrees to std vector of integer values
	std::vector<int> degrees = as<std::vector<int> >(pol_degrees);

	// Initiale degrees related variables
	int degrees_size = degrees.size();
	std::vector<int> degrees_products(degrees_size);

	// Calculate number of coefficients and degrees products
	
	int coefficients_size = 1;

	for (int i = 0; i < degrees_size; i++)
	{
		coefficients_size *= (degrees[i] + 1); //+1 because degrees order starts from zero
		degrees_products[i] = 1;

		for (int j = i + 1; j < degrees_size; j++)
		{
			degrees_products[i] *= (degrees[j] + 1);
		}
	}

	// Assign vector index to coefficients
	std::vector<std::vector<int>> ind_pattern_full = std::vector<std::vector<int>>(degrees_size);
	std::vector<std::vector<int>> coefficients_ind(coefficients_size);
	NumericMatrix coefficients_vec(degrees_size, coefficients_size);

	for (int i = 0; i < degrees_size; i++)
	{
		// Calculate pattern for i-th variable
		std::vector<int> ind_pattern = std::vector<int>(degrees_products[i] * (degrees[i] + 1));
		int counter = 0;

		for (int j = 0; j <= degrees[i]; j++)
		{
			for (int k = 0; k < degrees_products[i]; k++)
			{
				ind_pattern[counter] = j;
				counter++;
			}
		}

		int ind_pattern_times = coefficients_size / ind_pattern.size(); //times pattern repeats
		ind_pattern_full[i].reserve(coefficients_size); //preallocate memorry to increase insertation speed

		for (int j = 0; j < ind_pattern_times; j++)
		{
			// Repeat pattern untill the end of the pattern matrix row
			ind_pattern_full[i].insert(ind_pattern_full[i].end(), 
                                 ind_pattern.begin(), ind_pattern.end());
		}

		// Pattern defined for rows while coefficients indexes are located in columns
		// of pattern matrix. Lets copy values from patterns to coefficients indexes
		for (int j = 0; j < coefficients_size; j++)
		{
			coefficients_vec(i, j) = (double)(ind_pattern_full[i][j]);
		}
	}

	return(coefficients_vec);
}

//' Print polynomial given it's degrees and coefficients
//' @description This function prints polynomial given it's degrees and coefficients.
//' @template pol_degrees_Template
//' @template pol_coefficients_Template
//' @details Function automatically removes polynomial elements which coefficient are zero
//' and variables which power is zero. Output may contain long coefficients representation as they
//' are not rounded.
//' @return This function returns the string which contains polynomial symbolic representation.
//' @template printPol_examples_Template
//' @export
// [[Rcpp::export]]
Rcpp::String printPolynomial(NumericVector pol_degrees, NumericVector pol_coefficients)
{
	// Load R environments
	Rcpp::Environment base_env("package:base");
	Rcpp::Function paste0 = base_env["paste0"];

	// Get polynomial matrix from polynomialIndex function
	NumericVector pol_ind_mat = polynomialIndex(pol_degrees);

	// Create dimensions related variables
	int pol_coefficients_n = pol_coefficients.size();
	int pol_degrees_n = pol_degrees.size();

	// Initialize variable to store the polynomial symbolic representation
	std::string pol_string = "";

	// Iteratite throught polynomial coefficients and variables
	for (int i = 0; i < pol_coefficients_n; i++)
	{
		if ((pol_coefficients[i] != 0) | (i == 0))
		{
			if ((pol_coefficients[i] != 1) | (i == 0))
			{
				String pol_string_R = paste0((double)pol_coefficients[i]);
				pol_string += pol_string_R;
			}
			for (int j = 0; j < pol_degrees_n; j++)
			{
				int pol_pow = pol_ind_mat(j, i);
				if (pol_pow != 0)
				{
					pol_string += "x" + std::to_string(j + 1);
					if (pol_pow != 1)
					{
						pol_string += "^" + std::to_string(pol_pow);
					}
				}
			}
		}
			if (i < (pol_coefficients_n - 1))
			{
				if (pol_coefficients[i + 1] > 0)
				{
					pol_string += " + ";
				}
				if (pol_coefficients[i + 1] < 0)
				{
					pol_coefficients[i + 1] = -pol_coefficients[i + 1];
					pol_string += " - ";
				}
			}
	}

	return(pol_string);
}
