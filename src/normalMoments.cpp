#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends(RcppArmadillo)]]

//' Calculate k-th order moment of normal distribution
//' @description This function iteratively calculates k-th order moment of normal distribution.
//' @param k non-negative integer moment order.
//' @param mean numeric expected value.
//' @param sd positive numeric standard deviation.
//' @param return_all_moments logical; if \code{TRUE}, function returns the matrix (1 row, k+1 columns)
//' of moments of normaly distributed random variable with mean = \code{mean}
//' and standard deviation = \code{sd}. Note that i-th column value corresponds to the (i-1)-th moment.
//' @template is_validation_Template
//' @details This function estimates \code{k}-th order moment of
//' normal distribution which mean equals to \code{mean} and standard deviation equals to \code{sd}.\cr
//' @template k_integer_template
//' @return This function returns \code{k}-th order moment of
//' normal distribution which mean equals to \code{mean} and standard deviation is \code{sd}.
//' If \code{return_all_moments} is \code{TRUE} then see this argument description above for
//' output details.
//' @examples
//' ##Calculate 5-th order moment of normal random variable which
//' ##mean equals to 3 and standard deviation is 5.
//'
//' #5-th moment
//' normalMoment(k = 5, mean = 3, sd = 5)
//' 
//' #(0-5)-th moments
//' normalMoment(k = 5, mean = 3, sd = 5, return_all_moments = TRUE)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix normalMoment(int k = 0, 
						   double mean = 0, double sd = 1,
						   bool return_all_moments = false, 
						   bool is_validation = true)
{
	//Validation
	//------------------------------------------------------------
	if (is_validation)
	{
		if (k < 0)
		{
			stop("parameter k should be non-negative integer");
		}
		if (sd <= 0)
		{
			stop("parameter sd should be positive integer");
		}
	}
	//------------------------------------------------------------

	//Initialize matrix to store the moments
	NumericMatrix moments(1, k + 1);

	//If order is 0 just return 1
	if (k == 0)
	{
		moments[0] = 1;
		return(moments);
	}

	std::fill(moments.begin(), moments.end(), 1);

	moments[1] = mean;

	//If moment is 1 it's equals to mean
	if (k == 1)
	{
		if (!return_all_moments)
		{
			NumericMatrix moment = NumericMatrix(1, 1);
			moment(0, 0) = moments[k];
			return(moment);
		}
		return(moments);
	}

	//Iteratively calculate other moments

	double sd_squared = pow(sd, 2);

	for (int i = 2; i <= k; i++)
	{
		moments[i] = (i - 1) * sd_squared * moments[i - 2] + mean * moments[i - 1];
	}

	//Return depends on return_all_moments value

	if (!return_all_moments)
	{
		NumericMatrix moment = NumericMatrix(1, 1);
		moment(0, 0) = moments[k];
		return(moment);
	}

	return(moments);
}

//' Calculate k-th order moment of truncated normal distribution
//' @description This function iteratively calculates k-th order moment of truncated normal distribution.
//' @param k non-negative integer moment order.
//' @param x_lower numeric vector of lower trancation points.
//' @param x_upper numeric vector of upper trancation points.
//' @param mean numeric expected value.
//' @param sd positive numeric standard deviation.
//' @template pdf_lower_Template
//' @template cdf_lower_Template
//' @template pdf_upper_Template
//' @template cdf_upper_Template
//' @template cdf_difference_Template
//' @template is_validation_Template
//' @param return_all_moments logical; if \code{TRUE}, function returns the matrix of
//' moments of normaly distributed random variable with mean = \code{mean}
//' and standard deviation = \code{sd} under lower and upper truncation points
//' \code{x_lower} and \code{x_upper} correspondingly. Note that element in i-th row and
//' j-th column of this matrix corresponds to the i-th observation (j-1)-th
//' order moment.
//' @details This function estimates \code{k}-th order moment of
//' normal distribution which mean equals to \code{mean} and standard deviation equals to \code{sd} truncated
//' at points given by \code{x_lower} and \code{x_upper}. Note that the function is vectorized so you can provide
//' \code{x_lower} and \code{x_upper} as vectors of equal size. If vectors values for \code{x_lower} and \code{x_upper} are not
//' provided then their default values will be set to (-9999999999999.1) and (9999999999999.1) correspondingly.
//' @template k_integer_Template
//' @template pdf_cdf_precalculated_Template
//' @return This function returns vector of k-th order moments for normaly distributed random variable
//' with mean = \code{mean} and standard deviation = \code{sd} under x_lower
//' and x_upper truncation points \code{x_lower} and \code{x_upper} correspondingly.
//' If \code{return_all_moments} is \code{TRUE} then see this argument description above for
//' output details.
//' @examples
//' ##Calculate 5-th order moment of three truncated normal random variables (x1,x2,x3) 
//' ##which mean is 5 and standard deviation is 3. 
//' ##These random variables truncation points are given as follows:-1<x1<1, 0<x2<2, 1<x3<3.
//' k <- 3
//' x_lower <- c(-1, 0, 1)
//' x_upper <- c(1, 2 ,3)
//' mean <- 3
//' sd <- 5
//' 
//' #get the moments
//' truncatedNormalMoment(k, x_lower, x_upper, mean, sd)
//'
//' #get matrix of (0-5)-th moments (columns) for each variable (rows)
//' truncatedNormalMoment(k, x_lower, x_upper, mean, sd, return_all_moments = TRUE)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix truncatedNormalMoment(int k = 1,
	NumericVector x_lower = NumericVector(0),
	NumericVector x_upper = NumericVector(0),
	double mean = 0, double sd = 1,
	NumericVector pdf_lower = NumericVector(0),
	NumericVector cdf_lower = NumericVector(0),
	NumericVector pdf_upper = NumericVector(0),
	NumericVector cdf_upper = NumericVector(0),
	NumericVector cdf_difference = NumericVector(0),
	bool return_all_moments = false, 
	bool is_validation = true) 
{
	//Assign default truncation values

	if (x_lower.size() == 0)
	{
		x_lower = NumericVector::create(-9999999999999.1);
	}

	if (x_upper.size() == 0)
	{
		x_lower = NumericVector::create(9999999999999.1);
	}

	//Get number of observations
	int n = x_lower.size();

	//Validation
	//------------------------------------------------------------
	if (is_validation)
	{
		if (k < 0)
		{
			stop("parameter k should be non-negative integer");
		}
		if (sd <= 0)
		{
			stop("parameter sd should be positive integer");
		}
		if ((x_upper.size() != n) & (x_upper[0] != 9999999999999.1))
		{
			stop("vectors x_lower and x_upper should have the same length");
		}
		if ((x_lower.size() != n) & (x_lower[0] != -9999999999999.1))
		{
			stop("vectors x_lower and x_upper should have the same length");
		}
		if ((pdf_lower.size() != n) & (pdf_lower.size() != 0))
		{
			stop("vectors x_lower and pdf_lower should have the same length");
		}
		if ((cdf_lower.size() != n) & (cdf_lower.size() != 0))
		{
			stop("vectors x_lower and cdf_lower should have the same length");
		}
		if ((pdf_upper.size() != n) & (pdf_upper.size() != 0))
		{
			stop("vectors x_lower and pdf_upper should have the same length");
		}
		if ((cdf_upper.size() != n) & (cdf_upper.size() != 0))
		{
			stop("vectors x_lower and cdf_upper should have the same length");
		}
		if ((cdf_difference.size() != n) & (cdf_difference.size() != 0))
		{
			stop("vectors x_lower and cdf_difference should have the same length");
		}
	}
	//------------------------------------------------------------

	//Initialize matrix to store the moments
	NumericMatrix tr_moments(n, k + 1);

	std::fill(tr_moments.begin(), tr_moments.end(), 1);

	//If order is 0 just return 1
	if (k == 0)
	{
		return(tr_moments);
	}

	//PDF calculation (if not provided)
	if (pdf_lower.size() == 0)
	{
		pdf_lower = dnorm(x_lower, mean, sd);
	}
	if (pdf_upper.size() == 0)
	{
		pdf_upper = dnorm(x_upper, mean, sd);
	}

	//CDF calculation (if not provided)
	if (cdf_difference.size() == 0)
	{
		if (cdf_lower.size() == 0)
		{
			cdf_lower = pnorm(x_lower, mean, sd);
		}
		if (cdf_upper.size() == 0)
		{
			cdf_upper = pnorm(x_upper, mean, sd);
		}
		cdf_difference = cdf_upper - cdf_lower;
	}

	double sd_squared = pow(sd, 2);

	//The first moment
	tr_moments(_, 1) = (mean - sd_squared * ((pdf_upper - pdf_lower) / cdf_difference));

	if (k == 1)
	{
		return(tr_moments);
	}

	//Set infinity to zero in order to nullify power*pdf
	LogicalVector lower_cond = is_infinite(x_lower);
	LogicalVector upper_cond = is_infinite(x_upper);
	x_lower[lower_cond] = 0;
	x_upper[upper_cond] = 0;

	//Iteratively calculate other moments
	for (int i = 2; i <= k; i++)
	{
		tr_moments(_, i) = (i - 1) * sd_squared * tr_moments(_, i - 2) +
			mean * tr_moments(_, i - 1) - sd_squared * ((pow(x_upper, i - 1) * pdf_upper -
				pow(x_lower, i - 1) * pdf_lower) / cdf_difference);
	}

	//If return_all_moments is TRUE then return matrix of all moments from 0 to k
	if (return_all_moments)
	{
		return(tr_moments);
	}

	//If return_all_moments is FALSE then return k-th moment only
	NumericMatrix tr_moments_new(n, 1);
	tr_moments_new(_, 0) = tr_moments(_, k);

	return(tr_moments_new);
}
