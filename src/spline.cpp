#include "spline.h"
#include "normalMoments.h"
#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace RcppParallel;

//' @export
// [[Rcpp::export]]
List bsplineMult(List b,
                 double t1, 
                 double t2, 
                 bool is_left = true)
{
  // Prevent global environment changes
  b = clone(b);

  // Get the matrix of spline
  NumericMatrix b_m = b["m"];
  
  // Get the number of rows and columns of spline matrix
  int n_row = b_m.nrow();
  int n_col = b_m.ncol();

  // If the knot has multiplicity greater then one
  if(t1 == t2)
  {
    b_m = NumericMatrix(n_row, n_col + 1);
    b["m"] = b_m;
    
    return(b);
  }
  
  // Create matrix for the new spline
  NumericMatrix b_new = NumericMatrix(n_row, n_col + 1);
  
  // Estimate left or right part of the new spline
  if (is_left)
  {
    for (int i = 0; i < n_col; i++)
    {
      b_new(_, i + 1) = b_m(_, i);
      b_new(_, i) = b_new(_, i) - b_m(_, i) * t1;
    }
  } else {
    for (int i = 0; i < n_col; i++)
    {
      b_new(_, i + 1) = -b_m(_, i);
      b_new(_, i) = b_m(_, i) * t2 + b_new(_, i);
    }
  }
  b_new = b_new / (t2 - t1);
  
  // Insert a matrix into the spline
  b["m"] = b_new;
  
  return(b);
}

//' @export
// [[Rcpp::export]]
List bsplineMerge(List b_left, List b_right)
{
  // Get splines characteristics
  NumericVector knots = b_left["knots"];
  NumericVector ind = b_left["ind"];
  int i = ind[0];
  int d = ind[1] + 1;
  
  // Estimate left and right parts of the new spline
  List b1 = bsplineMult(b_left,
                        knots[i], knots[i + d],
                        true);
  List b2 = bsplineMult(b_right,
                        knots[i + 1], knots[i + d + 1],
                        false);

  // Combine the parts
  NumericMatrix b1_m = b1["m"];
  NumericMatrix b2_m = b2["m"];
  
  int n_col = b1_m.ncol();
  int n_row = b1_m.nrow();
  
  NumericMatrix b_m = NumericMatrix(n_row, n_col);
  for(int i = 0; i < n_col; i++)
  {
    b_m(_, i) = b1_m(_, i) + b2_m(_, i);
  }
  
  // Construct a result
  List b = clone(b1);
  b["m"] = b_m;
  b["ind"] = NumericVector::create(i, d);
  
  return(b);
}

//' @export
// [[Rcpp::export]]
List bsplineNames(List b)
{
  b = clone(b);
  NumericMatrix b_m = b["m"];
  NumericVector t = b["knots"];
  
  int n_row = b_m.nrow();
  int n_col = b_m.ncol();
  
  CharacterVector col_names = CharacterVector(n_col);
  for(int i = 0; i < n_col; i++)
  {
    col_names[i] = "x^" + std::to_string(i);
  }
  
  CharacterVector row_names = CharacterVector(n_row);
  for(int i = 0; i < (n_row - 1); i++)
  {
    row_names[i] = "[" + std::to_string(t[i]) + ", " + 
                         std::to_string(t[i + 1]) + ")";
  }
  row_names[n_row - 1] = "[" + std::to_string(t[n_row - 1]) + ", " + 
                               std::to_string(t[n_row]) + "]";
  
  rownames(b_m) = row_names;
  colnames(b_m) = col_names;
  
  return(b);
}

//' B-splines generation, estimation and combination
//' @name bspline
//' @description Function \code{\link[hpa]{bsplineGenerate}} generates a list
//' of all basis splines with appropriate \code{knots} vector and \code{degree}.
//' Function \code{\link[hpa]{bsplineComb}} allows to get linear combinations
//' of these b-splines with particular \code{weights}. 
//' Function \code{\link[hpa]{bsplineEstimate}} estimates the spline at
//' points \code{x}. The structure of this spline should be provided via
//' \code{m} and \code{knots} arguments.
//' @details In contrast to \code{\link[splines]{bs}} function 
//' \code{\link[hpa]{bsplineGenerate}} generates a splines basis in a form
//' of a list containing information concerning these b-splines structure.
//' In order to evaluate one of these b-splines at particular points
//' \code{\link[hpa]{bsplineEstimate}} function should be applied.
//' @return Function \code{\link[hpa]{bsplineGenerate}} returns a list. Each
//' element of this list is a list containing the following
//' information concerning b-spline structure:
//' \itemize{
//' \item \code{knots} - knots vector of the b-spline. 
//' \item \code{m} - matrix representing polynomial coefficients for each
//' interval of the spline in the same manner as for \code{m} argument
//' (see this argument description above).
//' \item \code{ind} - index of the b-spline.}
//' Function \code{bsplineComb} returns a list with the following arguments:
//' \itemize{
//' \item \code{knots} - knots vector of the \code{splines}. 
//' \item \code{m} - linear combination of the \code{splines} matrices; 
//' coefficients of this linear combination are given 
//' via \code{weights} argument.}
//' @return Function \code{\link[hpa]{bsplineGenerate}} returns a numeric
//' vector of values being calculated at points \code{x} via splines with 
//' \code{knots} vector and matrix \code{m}.
//' @template knots_Template
//' @template degree_Template
//' @template m_Template
//' @param is_names logical; if TRUE (default) then rows and columns of the
//' spline matrices will have a names. Set it to FALSE in order to get notable 
//' speed boost.
//' @template bsplines_examples_Template
//' @export
// [[Rcpp::export]]
List bsplineGenerate(NumericVector knots, 
                     int degree, 
                     bool is_names = true)
{
  // Get the number of knots
  int n_knots = knots.size();
  
  // Validate the input
  if (degree <= 0)
  {
    stop("degree should be positive integer");
  }
  
  for (int i = 0; i < (n_knots - 1); i++)
  {
    if(knots[i] > knots[i + 1])
    {
      stop("knots should be a vector of nondecreasing values");
    }
  }
  
  if (n_knots <= (degree + 1))
  {
    stop("The number of knots should be greater than (degree + 1)");
  }
  
  // Initial matrix for splines
  NumericMatrix init_mat = NumericMatrix(n_knots - 1, 1);
  init_mat(0, 0) = 1;
  
  // Fill initial splines
  List b(n_knots - 1);

  b[0] = List::create(Named("knots") = knots,
                      Named("ind") = NumericVector::create(0, 0),
                      Named("m") = clone(init_mat));

  for(int i = 1; i < (n_knots - 1); i++)
  {
    init_mat(i - 1, 0) = 0;
    init_mat(i, 0) = 1;

    b[i] = List::create(Named("knots") = knots,
                        Named("ind") = NumericVector::create(i, 0),
                        Named("m") = clone(init_mat));
  }

  // Calculate necessary splines
  for(int j = 1; j <= degree; j++)
  {
    List b_tmp(n_knots - 1 - j);
    for(int i = 0; i < (n_knots - 1 - j); i++)
    {
      b_tmp[i] = bsplineMerge(b[i], b[i + 1]);
    }
    b = b_tmp;
  }
  
  // Assign names
  if(is_names)
  {
    for(int i = 0; i < (n_knots - 1 - degree); i++)
    {
      b[i] = bsplineNames(b[i]);
    }
  }
    
  return(b);
}

//' @name bspline
//' @param x numeric vector representing the points at which the 
//' spline should be estimated.
//' @export
// [[Rcpp::export]]
NumericVector bsplineEstimate(NumericVector x,
                              NumericMatrix m, 
                              NumericVector knots)
{
  // Sore dimensions information
  int n = x.size();
  int m_col = m.ncol();
  int m_row = m.nrow();
  
  // Get degrees of the spline
  int degree = m_col - 1;
  
  // Estimate the matrix of
  // argument powers
  arma::mat m_arma = as<arma::mat>(m);
  arma::mat x_pow = arma::mat(n, m_col);
  arma::vec x_arma = as<arma::vec>(x);
  x_pow.col(0).fill(1);
  for(int i = 1; i <= degree; i++)
  {
    x_pow.col(i) = pow(x_arma, i);
  }

  // Calculate the final result
  arma::vec val = arma::vec(n);
  for (int i = 0; i < n; i++)
  {
    for(int j = 0; j < m_row; j++)
    {
      if ((knots[j] <= x[i]) & (x[i] <= knots[j + 1]))
      {
        val(i) = dot(m_arma.row(j), x_pow.row(i));
        break;
      }
    }
  }
  
  return(wrap(val));
}

//' @name bspline
//' @param splines list being returned by the 
//' \code{\link[hpa]{bsplineGenerate}} function or a manually constructed
//' list with b-splines knots and matrices entries.
//' @param weights numeric vector of the same length as \code{splines}.
//' @export
// [[Rcpp::export]]
List bsplineComb(List splines, 
                 NumericVector weights)
{
  // Get some variables
  List spline_new = splines[0];
  NumericMatrix m = spline_new["m"];
  m = clone(m);
  int n_splines = splines.size();
  int m_col = m.ncol();

  // Estimate new spline matrix
  m = m * weights[0];
  for(int i = 1; i < n_splines; i++)
  {
    List spline_tmp = splines[i];
    NumericMatrix m_tmp = spline_tmp["m"];
    for (int j = 0; j < m_col; j++)
    {
      m(_, j) = m(_, j) + m_tmp(_, j) * weights[i];
    }
  }
  
  // Organize and return the results
  List b = List::create(Named("knots") = spline_new["knots"],
                        Named("m") = m);
  
  return(b);
}

//' Probabilities and Moments Hermite Spline Approximation
//' @name hsaDist
//' @param x numeric vector of values for which the function should 
//' be estimated.
//' @template m_Template
//' @template knots_Template
//' @param mean expected value of a normal distribution.
//' @param sd standard deviation of a normal distribution.
//' @template log_Template
//' @description The set of functions similar to \code{\link[hpa]{dhpa}}-like
//' functions. The difference is that instead of polynomial these functions
//' utilize spline.
//' @details In contrast to \code{\link[hpa]{dhpa}}-like functions these
//' functions may deal with univariate distributions only. In future this
//' functions will be generalized to work with multivariate distributions.
//' The main idea of these functions is to use squared spline instead of squared 
//' polynomial in order to provide greater numeric stability and approximation 
//' accuracy. To provide spline parameters please use \code{m} and \code{knots}
//' arguments (i.e. instead of \code{pol_degrees} and \code{pol_coefficients}
//' arguments that where used to specify the polynomial
//' for \code{\link[hpa]{dhpa}}-like functions).
//' @return Function \code{\link[hpa]{dhsa}} returns vector of probabilities
//' of the same length as \code{x}. Function \code{\link[hpa]{ehsa}} 
//' returns moment value.
//' @seealso \code{\link[hpa]{dhpa}}, \code{\link[hpa]{bsplineGenerate}}
//' @template dhsa_examples_Template
//' @export
// [[Rcpp::export]]
NumericVector dhsa(NumericVector x,
                   NumericMatrix m,
                   NumericVector knots,
                   double mean = 0,
                   double sd = 1,
                   bool log = false)
{
  NumericVector val1 = bsplineEstimate(x, m, knots);

  if (log)
  {
    val1 = 2 * Rcpp::log(Rcpp::abs(val1)) + dnorm(x, mean, sd, true);
  } else {
    val1 = pow(val1, 2) * dnorm(x, mean, sd);
  }

  int degree = m.ncol() - 1;
  int n_knots = knots.size();
  
  Range r1 = Range(0, n_knots - 2);
  Range r2 = Range(1, n_knots - 1);
  
  NumericVector cdf_knots = pnorm(knots, mean, sd);
  NumericVector cdf_diff = cdf_knots[r2] - cdf_knots[r1];
  
  NumericMatrix moments = truncatedNormalMoment(2 * degree,  
                                                knots[r1], 
                                                knots[r2], 
                                                mean, sd, 
                                                NumericVector(0), cdf_knots[r1],
                                                NumericVector(0), cdf_knots[r2],
                                                cdf_diff,
                                                true, false, false, "NO");

  double val2 = 0;
  for (int t = 0; t < (n_knots - 1); t++)
  {
    if (cdf_diff[t] != 0)
    {
      for (int i = 0; i <= degree; i++)
      {
        for (int j = 0; j <= degree; j++)
        {
          val2 = val2 + m(t, i) * m(t, j) *
                        moments(t, i + j) * cdf_diff[t];
        }
      }
    }
  }

  NumericVector val;
  if (log)
  {
    val = val1 - std::log(val2);  
  } else {
    val = val1 / val2;
  }
  
  return(val);
}

//' @name hsaDist
//' @param power non-negative integer representing the power of the 
//' expected value i.e. E(X ^ power) will be estimated.
//' @export
// [[Rcpp::export]]
double ehsa(NumericMatrix m,
            NumericVector knots,
            double mean = 0,
            double sd = 1,
            double power = 1)
{
  int degree = m.ncol() - 1;
  int n_knots = knots.size();
  
  Range r1 = Range(0, n_knots - 2);
  Range r2 = Range(1, n_knots - 1);
  
  NumericVector cdf_knots = pnorm(knots, mean, sd);
  NumericVector cdf_diff = cdf_knots[r2] - cdf_knots[r1];
  
  NumericMatrix moments = truncatedNormalMoment(2 * degree + power,  
                                                knots[r1], knots[r2], 
                                                mean, sd, 
                                                NumericVector(0), cdf_knots[r1],
                                                NumericVector(0), cdf_knots[r2],
                                                cdf_diff,
                                                true, false, false, "NO");
  
  double val1 = 0;
  double val2 = 0;
  for (int t = 0; t < (n_knots - 1); t++)
  {
    if (cdf_diff[t] != 0)
    {
      for (int i = 0; i <= degree; i++)
      {
        for (int j = 0; j <= degree; j++)
        {
          double m_prod = m(t, i) * m(t, j);
          val1 = val1 + m_prod * cdf_diff[t] * moments(t, i + j + power);
          val2 = val2 + m_prod * cdf_diff[t] * moments(t, i + j);
        }
      }
    }
  }
  
  double val = val1 / val2;
  
  return(val);
}
