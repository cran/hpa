#ifndef hpa_spline_H
#define hpa_spline_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;


List bsplineMult(List b, 
                 double t1, 
                 double t2, 
                 bool is_left);

List bsplineMerge(List b_left,
                  List b_right);

List bsplineNames(List b);

List bsplineGenerate(NumericVector knots, 
                     int degree,
                     bool is_names);

NumericVector bsplineEstimate(NumericVector x, 
                              NumericVector knots,
                              NumericMatrix m);

List bsplineComb(List splines, 
                 NumericVector weights);

NumericVector dhsa(NumericVector x, 
                   NumericMatrix m, 
                   NumericVector knots, 
                   double mean, 
                   double sd, 
                   bool log);

double ehsa(NumericMatrix m, 
            NumericVector knots, 
            double mean, 
            double sd,
            double power);

#endif
