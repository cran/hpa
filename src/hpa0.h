#ifndef hpa_hpa0_H
#define hpa_hpa0_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;

List dhpa0(
    const arma::vec x,
    const arma::vec pc,
    double mean,
    double sd,
    bool is_parallel,
    bool log,
    bool is_validatione,
    bool is_grad);

List phpa0(
    const arma::vec x,
    const arma::vec pc,
    double mean,
    double sd,
    bool is_parallel,
    bool log,
    bool is_validatione,
    bool is_grad);

#endif
