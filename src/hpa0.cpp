#define ARMA_DONT_USE_OPENMP
#include <RcppArmadillo.h>
#include "normalMoments.h"

using namespace Rcpp;
using namespace RcppArmadillo;

// [[Rcpp::interfaces(r, cpp)]]

//' Fast pdf and cdf for standardized univariate PGN distribution
//' @description This function uses fast algorithms to calculate densities
//' and probabilities (along with their derivatives) related to standardized 
//' PGN distribution.
//' @name hpaDist0
//' @param x numeric vector of functions arguments.
//' @param pc polynomial coefficients without the first term.
//' @param mean expected value (mean) of the distribution.
//' @param sd standard deviation of the distribution.
//' @param is_parallel logical; if TRUE then multiple cores will be used for 
//' some calculations. Currently unavailable.
//' @template log_Template
//' @template is_validation_Template
//' @param is_grad logical; if \code{TRUE} (default) then function returns 
//' gradients respect to \code{x} and \code{pc}.
//' @details Functions \code{\link[hpa]{dhpa0}} and 
//' \code{\link[hpa]{phpa0}} are similar to \code{\link[hpa]{dhpa}} and
//' \code{\link[hpa]{phpa}} correspondingly. However there are two key
//' differences. First, \code{\link[hpa]{dhpa0}} and \code{\link[hpa]{phpa0}}
//' are deal with univariate PGN distribution only. Second, this distribution
//' is standardized to zero mean and unit variances. Moreover \code{pc} is 
//' similar to \code{pol_coefficients} argument of \code{\link[hpa]{dhpa}} but
//' without the first component i.e. \code{pc=pol_coefficients[-1]}. Also
//' \code{mean} and \code{sd} are not the arguments of the normal density
//' but actual mean and standard deviation of the resulting distribution. So
//' if these arguments are different from \code{0} and \code{1} correspondingly
//' then standardized PGN distribution will be linearly transformed to have
//' mean \code{mean} and standard deviation \code{sd}.
//' @return Both functions return a list.
//' Function \code{\link[hpa]{dhpa0}} returns a list with element named
//' \code{"den"} that is a numeric vector of density values. 
//' Function \code{\link[hpa]{phpa0}} returns a list with element named
//' \code{"prob"} that is a numeric vector of probabilities. 
//' 
//' If \code{is_grad = TRUE} then elements \code{"grad_x"} and \code{"grad_pc"}
//' will be add to the list containing gradients respect to input argument
//' \code{x} and parameters \code{pc} correspondingly. If \code{log = TRUE} then
//' additional elements will be add to the list containing density, probability
//' and gradient values for logarithms of corresponding functions. These
//' elements will be named as \code{"grad_x_log"}, \code{"grad_pc_log"},
//' \code{"grad_prob_log"} and \code{"grad_den_log"}.
//' @examples
//' # Calculate density and probability of standartized PGN
//' # distribution
//'   # distribution parameters
//' pc <- c(0.5, -0.2)
//'   # function arguments
//' x <- c(-0.3, 0.8, 1.5)
//'   # probability density function
//' dhpa0(x, pc)
//'   # cumulative distribution function
//' phpa0(x, pc)
//' 
//' # Additionally calculate gradients respect to arguments
//' # and parameters of the PGN distribution
//' dhpa0(x, pc, is_grad = TRUE)
//' phpa0(x, pc, is_grad = TRUE)
//' 
//' # Let's denote by X standardized PGN random variable and repeat
//' # calculations for 2 * X + 1
//' dhpa0(x, pc, is_grad = TRUE, mean = 1, sd = 2)
//' phpa0(x, pc, is_grad = TRUE, mean = 1, sd = 2)
//' @export
// [[Rcpp::export(rng = false)]]
List dhpa0(
    const arma::vec x,
    const arma::vec pc,
    double mean = 0,
    double sd = 1,
    bool is_parallel = false,
    bool log = false,
    bool is_validation = true,
    bool is_grad = false)
{
  // Validation if need
  if (is_validation)
  {
    if (sd <= 0)
    {
      stop("Parameter 'sd' should be positive.");
    }
    
    if (pc.size() > 15)
    {
      stop("Currently length of 'pc' should not be greater than 15.");
    }
  }
  
  // Create output list
  List return_list;
  
  // Some dimensions related constants
  const int n = x.size();
  const int K = pc.size();
  
  // Deal with infinite values if need
  if (x.has_inf())
  {
    // Find infinite and finite values indexes
    arma::uvec finite_ind = arma::find_finite(x);
    arma::uvec infinite_neg_ind = arma::find(x == (-arma::datum::inf));
    arma::uvec infinite_ind = arma::find(x == arma::datum::inf);
    arma::vec x_finite = x.elem(finite_ind);
    
    // Calculate density for infinite values
    arma::vec den_new = arma::vec(n);
    
    // Take logarithm of infinite probabilities if need
    arma::vec den_log_new = arma::vec(n, arma::fill::value(-arma::datum::inf));
    
    // Prepare zero derivatives for infinite arguments
    arma::vec grad_x_new = arma::vec(n);
    arma::vec grad_x_log_new = arma::vec(n);
    arma::mat grad_pc_new = arma::mat(n, K);
    arma::mat grad_pc_log_new = arma::mat(n, K);
    
    // Get estimates for finite values
    List dhpa0_val;
    if (finite_ind.size() > 0)
    {
      dhpa0_val = dhpa0(x_finite, pc, 
                        mean, sd, 
                        is_parallel, log, 
                        false, is_grad);
    }
    else
    {
      return_list["den"] = den_new;
      return_list["den_log"] = den_log_new;
      return_list["grad_x"] = grad_x_new;
      return_list["grad_x_log"] = grad_x_log_new;
      return_list["grad_pc"] = grad_pc_new;
      return_list["grad_pc_log"] = grad_pc_log_new;
      return(return_list);
    }
    
    // Combine estimates for finite and infinite values
    if (dhpa0_val.containsElementNamed("den"))
    {
      arma::vec den_new0 = dhpa0_val["den"];
      den_new.elem(finite_ind) = den_new0;
      return_list["den"] = den_new;
    }
    
    if (dhpa0_val.containsElementNamed("den_log"))
    {
      arma::vec den_log_new0 = dhpa0_val["den_log"];
      den_log_new.elem(finite_ind) = den_log_new0;
      return_list["den_log"] = den_log_new;
    }
    
    if (dhpa0_val.containsElementNamed("grad_x"))
    {
      arma::vec grad_x_new0 = dhpa0_val["grad_x"];
      grad_x_new.elem(finite_ind) = grad_x_new0;
      return_list["grad_x"] = grad_x_new;
    }
    
    if (dhpa0_val.containsElementNamed("grad_x_log"))
    {
      arma::vec grad_x_log_new0 = dhpa0_val["grad_x_log"];
      grad_x_log_new.elem(finite_ind) = grad_x_log_new0;
      return_list["grad_x_log"] = grad_x_log_new;
    }
    
    if (dhpa0_val.containsElementNamed("grad_pc"))
    {
      arma::mat grad_pc_new0 = dhpa0_val["grad_pc"];
      grad_pc_new.rows(finite_ind) = grad_pc_new0;
      return_list["grad_pc"] = grad_pc_new;
    }
    
    if (dhpa0_val.containsElementNamed("grad_pc_log"))
    {
      arma::mat grad_pc_log_new0 = dhpa0_val["grad_pc_log"];
      grad_pc_log_new.rows(finite_ind) = grad_pc_log_new0;
      return_list["grad_pc_log"] = grad_pc_log_new;
    }
    
    return(return_list);
  }
  
  // Polynomial coefficients with one
  arma::vec pc1(K + 1);
  pc1.at(0) = 1;
  pc1.subvec(1, K) = pc;
  
  // Denominator and moments nominator
  arma::vec moments = {1, 0, 1, 0, 3, 0, 15, 0, 105, 0, 945, 0, 10395, 0,
                       135135, 0, 2027025, 0, 34459425, 0, 654729075};
  double denom = 0;
  double e1_nom = 0;
  double e2_nom = 0;
  for (int i = 0; i <= K; i++)
  {
    for (int j = i; j <= K; j++)
    {
      int sum_i_j = i + j;
      double pc_prod = pc1.at(i) * pc1.at(j);
      // density and first moment
      if (sum_i_j % 2 == 0)
      {
        if (i != j)
        {
          denom = denom + 2 * pc_prod * moments.at(sum_i_j);
          e2_nom = e2_nom + 2 * pc_prod * moments.at(sum_i_j + 2);
        }
        else
        {
          denom = denom + pc_prod * moments.at(sum_i_j);
          e2_nom = e2_nom + pc_prod * moments.at(sum_i_j + 2);
        }
      }
      // second moment
      else
      {
        if (i != j)
        {
          e1_nom = e1_nom + 2 * pc_prod * moments.at(sum_i_j + 1);
        }
        else
        {
          e1_nom = e1_nom + pc_prod * moments.at(sum_i_j + 1);
        }
      }
    }
  }
  
  // Calculate moments
  double e1 = e1_nom / denom;
  double e2 = e2_nom / denom;
  double evar = e2 - pow(e1, 2);
  double esd = sqrt(evar);
  double sd_adj = esd / sd;
  
  // Nominator
  arma::vec x0 = x - mean;
  arma::vec x_adj = sd_adj * x0 + e1;
  arma::mat x_pow(n, K + 1);
  x_pow.col(0) = arma::ones(n);
  x_pow.col(1) = x_adj;
  for (int i = 2; i <= K; i++)
  {
    x_pow.col(i) = x_pow.col(i - 1) % x_adj;
  }
  arma::vec x_pc = x_pow * pc1;
  arma::vec den_nom = arma::pow(x_pc, 2);

  // Normal density
  arma::vec den_norm = arma::normpdf(x_adj);
  
  // Density value
  arma::vec den = (sd_adj / denom) * (den_norm % den_nom);
  return_list["den"] = den;
  
  // Take the log if need
  arma::vec den_log;
  if (log)
  {
    den_log = arma::log(den);
    return_list["den_log"] = den_log;
  }
  
  if (!is_grad)
  {
    return(return_list);
  }
  
  // Calculate derivative respect to argument of the function
  arma::vec den_d_x_unadj = (x_pow.cols(0, K - 1) * 
                             ((2 * pc) % arma::linspace<arma::vec>(1, K, K))) / 
                            x_pc - x_adj;
  arma::vec den_d_x = den_d_x_unadj * sd_adj;
  return_list["grad_x_log"] = den_d_x;
  
  if (!log)
  {
    return_list["grad_x"] = den_d_x % den;
  }
  
  // Calculate moments related parts
  arma::vec e0_d_pc(K);
  arma::vec e1_d_pc(K);
  arma::vec e2_d_pc(K);
  for (int i = 0; i < K; i++)
  {
    e0_d_pc.at(i) = arma::dot(pc1, moments.subvec(i + 1, i + 1 + K));
    e1_d_pc.at(i) = arma::dot(pc1, moments.subvec(i + 2, i + 2 + K));
    e1_d_pc.at(i) = e1_d_pc.at(i) * denom - e1_nom * e0_d_pc.at(i);
    e2_d_pc.at(i) = arma::dot(pc1, moments.subvec(i + 3, i + 3 + K));
    e2_d_pc.at(i) = e2_d_pc.at(i) * denom - e2_nom * e0_d_pc.at(i);
  }
  double denom_sqr = pow(denom, 2);
  e0_d_pc = e0_d_pc / denom;
  e1_d_pc = (2 * e1_d_pc) / denom_sqr;
  e2_d_pc = (2 * e2_d_pc) / denom_sqr;
  arma::vec var_d_pc = e2_d_pc - 2 * e1 * e1_d_pc;
  arma::vec sd_d_pc = var_d_pc / (2 * esd);
  arma::vec e_d_pc = sd_d_pc / sd;
  
  // Calculate derivative respect to log-density before
  // accounting for moments
  arma::mat den_d_pc = x_pow.cols(1, K).each_col() / x_pc;
  den_d_pc = 2 * (den_d_pc.each_row() - e0_d_pc.t());
  
  // Account for moments
  for (int i = 0; i < K; i++)
  {
    den_d_pc.col(i) = den_d_pc.col(i) + 
                      ((e_d_pc.at(i) * x0 + e1_d_pc.at(i)) % den_d_x_unadj) +
                      (sd_d_pc.at(i) / esd);
   }
  return_list["grad_pc_log"] = den_d_pc;
  
  // // Derivative respect to standard deviation
  // arma::vec grad_sd_log = den_d_x_unadj % ((-sd_adj / sd) * x0) - (1 / sd);
  // return_list["grad_sd_log"] = grad_sd_log;
  // 
  // // Derivative respect to mean
  // arma::vec grad_mean_log = den_d_x_unadj * (-sd_adj);
  // return_list["grad_mean_log"] = grad_mean_log;
  
  // If not logarithm
  if (!log)
  {
    return_list["grad_pc"] = den_d_pc.each_col() % den;
    // return_list["grad_sd"] = grad_sd_log % den;
    // return_list["grad_mean"] = grad_mean_log % den;
  }
  
  return(return_list);
}

//' @name hpaDist0
//' @export
// [[Rcpp::export(rng = false)]]
List phpa0(
    const arma::vec x,
    const arma::vec pc,
    double mean = 0,
    double sd = 1,
    bool is_parallel = false,
    bool log = false,
    bool is_validation = true,
    bool is_grad = false)
{
  // Validation if need
  if (is_validation)
  {
    if (sd <= 0)
    {
      stop("Parameter 'sd' should be positive.");
    }
    
    if (pc.size() > 15)
    {
      stop("Currently length of 'pc' should not be greater than 15.");
    }
  }
  
  // Create output list
  List return_list;
  
  // Some dimensions related constants
  const int n = x.size();
  const int K = pc.size();
  const int m = 2 * K + 1;
  
  // Deal with infinite values if need
  if (x.has_inf())
  {
    // Find infinite and finite values indexes
    arma::uvec finite_ind = arma::find_finite(x);
    arma::uvec infinite_neg_ind = arma::find(x == (-arma::datum::inf));
    arma::uvec infinite_ind = arma::find(x == arma::datum::inf);
    
    // Calculate probabilities for infinite values
    arma::vec x_finite = x.elem(finite_ind);
    arma::vec prob_new = arma::vec(n);
    prob_new.elem(infinite_ind).ones();
    
    // Take logarithm of infinite probabilities if need
    arma::vec prob_log_new = arma::vec(n);
    prob_log_new.elem(infinite_neg_ind).fill(-arma::datum::inf);
    
    // Prepare zero derivatives for infinite arguments
    arma::vec grad_x_new = arma::vec(n);
    arma::vec grad_x_log_new = arma::vec(n);
    arma::mat grad_pc_new = arma::mat(n, K);
    arma::mat grad_pc_log_new = arma::mat(n, K);
    
    // Get estimates for finite values
    List phpa0_val;
    if (finite_ind.size() > 0)
    {
      phpa0_val = phpa0(x_finite, pc, 
                        mean, sd, 
                        is_parallel, log, 
                        false, is_grad);
    }
    else
    {
      return_list["prob"] = prob_new;
      return_list["prob_log"] = prob_log_new;
      return_list["grad_x"] = grad_x_new;
      return_list["grad_x_log"] = grad_x_log_new;
      return_list["grad_pc"] = grad_pc_new;
      return_list["grad_pc_log"] = grad_pc_log_new;
      return(return_list);
    }
    
    // Combine estimates for finite and infinite values
    if (phpa0_val.containsElementNamed("prob"))
    {
      arma::vec prob_new0 = phpa0_val["prob"];
      prob_new.elem(finite_ind) = prob_new0;
      return_list["prob"] = prob_new;
    }
    
    if (phpa0_val.containsElementNamed("prob_log"))
    {
      arma::vec prob_log_new0 = phpa0_val["prob_log"];
      prob_log_new.elem(finite_ind) = prob_log_new0;
      return_list["prob_log"] = prob_log_new;
    }
    
    if (phpa0_val.containsElementNamed("grad_x"))
    {
      arma::vec grad_x_new0 = phpa0_val["grad_x"];
      grad_x_new.elem(finite_ind) = grad_x_new0;
      return_list["grad_x"] = grad_x_new;
    }
    
    if (phpa0_val.containsElementNamed("grad_x_log"))
    {
      arma::vec grad_x_log_new0 = phpa0_val["grad_x_log"];
      grad_x_log_new.elem(finite_ind) = grad_x_log_new0;
      return_list["grad_x_log"] = grad_x_log_new;
    }
    
    if (phpa0_val.containsElementNamed("grad_pc"))
    {
      arma::mat grad_pc_new0 = phpa0_val["grad_pc"];
      grad_pc_new.rows(finite_ind) = grad_pc_new0;
      return_list["grad_pc"] = grad_pc_new;
    }
    
    if (phpa0_val.containsElementNamed("grad_pc_log"))
    {
      arma::mat grad_pc_log_new0 = phpa0_val["grad_pc_log"];
      grad_pc_log_new.rows(finite_ind) = grad_pc_log_new0;
      return_list["grad_pc_log"] = grad_pc_log_new;
    }
    
    return(return_list);
  }
  
  // Polynomial coefficients with one
  arma::vec pc1(K + 1);
  pc1.at(0) = 1;
  pc1.subvec(1, K) = pc;
  
  // Denominator and moments nominator
  arma::vec moments = {1, 0, 1, 0, 3, 0, 15, 0, 105, 0, 945, 0, 10395, 0,
                       135135, 0, 2027025, 0, 34459425, 0, 654729075, 0,
                       13749310575, 0, 316234143225, 0, 7905853580625, 0,
                       213458046676875, 0, 6190283353629375, 0};
  double denom = 0;
  double e1_nom = 0;
  double e2_nom = 0;
  for (int i = 0; i <= K; i++)
  {
    for (int j = i; j <= K; j++)
    {
      int sum_i_j = i + j;
      double pc_prod = pc1.at(i) * pc1.at(j);
      // density and first moment
      if (sum_i_j % 2 == 0)
      {
        if (i != j)
        {
          denom = denom + 2 * (pc_prod * moments.at(sum_i_j));
          e2_nom = e2_nom + 2 * (pc_prod * moments.at(sum_i_j + 2));
        }
        else
        {
          denom = denom + pc_prod * moments.at(sum_i_j);
          e2_nom = e2_nom + pc_prod * moments.at(sum_i_j + 2);
        }
      }
      // second moment
      else
      {
        if (i != j)
        {
          e1_nom = e1_nom + 2 * pc_prod * moments.at(sum_i_j + 1);
        }
        else
        {
          e1_nom = e1_nom + pc_prod * moments.at(sum_i_j + 1);
        }
      }
    }
  }
  
  // Calculate moments
  double e1 = e1_nom / denom;
  double e2 = e2_nom / denom;
  double evar = e2 - pow(e1, 2);
  double esd = sqrt(evar);
  double sd_adj = esd / sd;
  
  // Adjust the argument
  arma::vec x0 = x - mean;
  arma::vec x_adj = sd_adj * x0 + e1;
  
  // Normal density and probability
  arma::vec den_norm = arma::normpdf(x_adj);
  arma::vec prob_norm = arma::normcdf(x_adj);
  arma::vec neg_ratio = -den_norm / prob_norm;

  // Powers
  arma::mat x_pow(n, m);
  x_pow.col(0) = arma::ones(n);
  x_pow.col(1) = x_adj;
  for (int i = 2; i < m; i++)
  {
    x_pow.col(i) = x_pow.col(i - 1) % x_adj;
  }
  
  // Truncated normal moments
  arma::mat tr_moments(n, m);
  tr_moments.col(0) = arma::ones(n);
  tr_moments.col(1) = neg_ratio;
  for (int i = 2; i < m; i++)
  {
    tr_moments.col(i) = (i - 1) * tr_moments.col(i - 2) + 
                        neg_ratio % x_pow.col(i - 1);
      
  }
  
  // Calculate the nominator
  arma::vec nom(n);
  for (int i = 0; i <= K; i++)
  {
    for (int j = i; j <= K; j++)
    {
      if (i == j)
      {
        nom = nom + pow(pc1.at(i), 2) * tr_moments.col(i + j);
      }
      else
      {
        nom = nom + (2 * pc1.at(i) * pc1.at(j)) * tr_moments.col(i + j);
      }
    }
  }
  
  // Aggregate the results into probability
  arma::vec prob = (prob_norm % nom) / denom;
  return_list["prob"] = prob;
  
  arma::vec prob_log;
  if (log | is_grad)
  {
    prob_log = arma::log(prob);
    return_list["prob_log"] = prob_log;
  }
  
  if (!is_grad)
  {
    return(return_list);
  }
  
  // Calculate density if need
  arma::vec x_pc = x_pow.cols(0, K) * pc1;
  arma::vec den_nom = arma::pow(x_pc, 2);
  
  // Density value
  arma::vec den = (sd_adj / denom) * (den_norm % den_nom);
  return_list["grad_x"] = den;
  
  arma::vec grad_x_log = den / prob;
  return_list["grad_x_log"] = grad_x_log;
  
  // Calculate moments related parts
  arma::vec e0_d_pc(K);
  arma::vec e1_d_pc(K);
  arma::vec e2_d_pc(K);
  arma::mat etr_d_pc(n, K);
  for (int i = 0; i < K; i++)
  {
    e0_d_pc.at(i) = arma::dot(pc1, moments.subvec(i + 1, i + 1 + K));
    e1_d_pc.at(i) = arma::dot(pc1, moments.subvec(i + 2, i + 2 + K));
    e1_d_pc.at(i) = e1_d_pc.at(i) * denom - e1_nom * e0_d_pc.at(i);
    e2_d_pc.at(i) = arma::dot(pc1, moments.subvec(i + 3, i + 3 + K));
    e2_d_pc.at(i) = e2_d_pc.at(i) * denom - e2_nom * e0_d_pc.at(i);
    etr_d_pc.col(i) = tr_moments.cols(i + 1, i + 1 + K) * pc1;
  }
  double denom_sqr_adj = 2 / pow(denom, 2);
  e0_d_pc = e0_d_pc / denom;
  e1_d_pc = e1_d_pc * denom_sqr_adj;
  e2_d_pc = e2_d_pc * denom_sqr_adj;
  arma::vec var_d_pc = e2_d_pc - 2 * e1 * e1_d_pc;
  arma::vec sd_d_pc = var_d_pc / (2 * esd);
  arma::vec e_d_pc = sd_d_pc / sd;
  
  // Calculate derivative respect to log-probability before
  // accounting for moments
  arma::mat prob_d_pc = etr_d_pc.each_col() / nom;
  prob_d_pc = 2 * (prob_d_pc.each_row() - e0_d_pc.t());

  // Account for moments
  arma::vec grad_x_log_adj = grad_x_log / sd_adj;
  for (int i = 0; i < K; i++)
  {
    prob_d_pc.col(i) = prob_d_pc.col(i) + 
                       ((e_d_pc.at(i) * x0 + e1_d_pc.at(i)) % grad_x_log_adj);
  }
  return_list["grad_pc_log"] = prob_d_pc;
  
  // // Derivative respect to standard deviation
  // arma::vec grad_sd_log = grad_x_log_adj % ((-sd_adj / sd) * x0);
  // return_list["grad_sd_log"] = grad_sd_log;
  // 
  // // Derivative respect to mean
  // arma::vec grad_mean_log = grad_x_log_adj * (-sd_adj);
  // return_list["grad_mean_log"] = grad_mean_log;
  
  // If not logarithm
  if (!log)
  {
    return_list["grad_pc"] = prob_d_pc.each_col() % prob;
    // return_list["grad_sd"] = grad_sd_log % prob;
    // return_list["grad_mean"] = grad_mean_log % prob;
  }
  
  return(return_list);
}
