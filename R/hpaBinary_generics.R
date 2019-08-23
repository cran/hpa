#' Predict method for hpaBinary
#' @param object Object of class "hpaBinary"
#' @template newdata_Template
#' @param is_prob logical; if TRUE (default) then function returns 
#' predicted probabilities. Otherwise latent variable
#' (single index) estimates will be returned.
#' @template elipsis_Template
#' @return This function returns predicted probabilities based on 
#' \code{\link[hpa]{hpaBinary}} estimation results.
predict.hpaBinary <- function (object, ..., 
                              newdata = NULL, 
                              is_prob = TRUE)
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  return(predict_hpaBinary(object, newdata, is_prob))
}
###
#' Summarizing hpaBinary Fits
#' @param object Object of class "hpaBinary"
#' @template elipsis_Template
#' @return This function returns the same list as \code{\link[hpa]{hpaBinary}} 
#' function changing it's class to "summary.hpaBinary".
summary.hpaBinary <- function (object, ...) 
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  return(summary_hpaBinary(object))
}
###
#' Summary for hpaBinary output
#' @param x Object of class "hpaML"
#' @template elipsis_Template
print.summary.hpaBinary <- function (x, ...) 
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  return(print_summary_hpaBinary(x))
}
###
#' Plot hpaBinary random errors approximated density
#' @param x Object of class "hpaBinary"
#' @template elipsis_Template
plot.hpaBinary <- function (x, ...) 
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  return(plot_hpaBinary(x))
}
###
#' Calculates AIC for "hpaBinary" object
#' @description This function calculates AIC for "hpaBinary" object
#' @param object Object of class "hpaBinary"
#' @template elipsis_Template
#' @template AIC_Template
AIC.hpaBinary <- function (object, ..., k = 2)
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  return(AIC_hpaBinary(object, k))
}
###
#' Calculates log-likelihood for "hpaBinary" object
#' @description This function calculates log-likelihood for "hpaBinary" object
#' @param object Object of class "hpaBinary"
#' @template elipsis_Template
logLik.hpaBinary <- function (object, ...)
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  return(logLik_hpaBinary(object))
}