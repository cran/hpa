#' Predict method for hpaML
#' @param object Object of class "hpaML"
#' @template newdata_Template
#' @template elipsis_Template
#' @return This function returns predictions based on 
#' \code{\link[hpa]{hpaML}} estimation results.
predict.hpaML <- function (object, ..., newdata = matrix(c(0)))
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  return(predict_hpaML(object, as.matrix(newdata)))
}
###
#' Summarizing hpaML Fits
#' @param object Object of class "hpaML"
#' @template elipsis_Template
#' @return This function returns the same list as \code{\link[hpa]{hpaML}} 
#' function changing it's class to "summary.hpaML".
summary.hpaML <- function (object, ...)
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  return(summary_hpaML(object))
}
###
#' Summary for hpaML output
#' @param x Object of class "hpaML"
#' @template elipsis_Template
print.summary.hpaML <- function (x, ...)
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  return(print_summary_hpaML(x))
}
###
#' Calculates AIC for "hpaML" object
#' @description This function calculates AIC for "hpaML" object
#' @param object Object of class "hpaML"
#' @template elipsis_Template
#' @template AIC_template
AIC.hpaML <- function (object, ..., k = 2)
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  return(AIC_hpaML(object))
}
###
#' Calculates log-likelihood for "hpaML" object
#' @description This function calculates log-likelihood for "hpaML" object
#' @param object Object of class "hpaML"
#' @template elipsis_Template
logLik.hpaML <- function (object, ...)
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  return(logLik_hpaML(object, ...))
}
