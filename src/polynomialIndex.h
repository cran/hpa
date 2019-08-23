#ifndef hpa_polynomialIndex_H
#define hpa_polynomialIndex_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;

	Rcpp::NumericMatrix polynomialIndex(Rcpp::NumericVector pol_degrees);

	Rcpp::String printPolynomial(Rcpp::NumericVector pol_degrees, 
								 Rcpp::NumericVector pol_coefficients);

#endif
