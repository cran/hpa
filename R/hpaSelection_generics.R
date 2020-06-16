#' Predict outcome and selection equation values from hpaSelection model
#' @description This function predicts outcome and selection equation values from hpaSelection model.
#' @param object Object of class "hpaSelection"
#' @param method string value indicating prediction method based on hermite polynomial approximation "HPA" or Newey method "Newey".
#' @template newdata_Template
#' @param is_cond logical; if \code{TRUE} (default) then conditional predictions will be estimated. Otherwise unconditional predictions will be returned.
#' @param is_outcome logical; if \code{TRUE} (default) then predictions for selection equation will be estimated using "HPA" method.
#' Otherwise selection equation predictions (probabilities) will be returned.
#' @template elipsis_Template
#' @details Note that Newey method can't predict conditional outcomes for zero selection equation value. Conditional probabilities for selection equation
#' could be estimated only when dependent variable from outcome equation is observable.
#' @return This function returns the list which structure depends on \code{method}, \code{is_probit} and \code{is_outcome} values.
predict.hpaSelection <- function (object,  ...,
                              newdata = NULL, method = "HPA", 
                              is_cond = TRUE, 
                              is_outcome = TRUE)
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  return(predict_hpaSelection(object, newdata, 
                       method, is_cond, 
                       is_outcome))
}
###
#' Summarizing hpaSelection Fits
#' @description This function summarizing hpaSelection Fits
#' @param object Object of class "hpaSelection"
#' @template elipsis_Template
#' @return This function returns the same list as \code{\link[hpa]{hpaSelection}} 
#' function changing it's class to "summary.hpaSelection".
summary.hpaSelection <- function (object, ...) 
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  return(summary_hpaSelection(object))
}
###
#' Summary for hpaSelection output
#' @param x Object of class "hpaSelection"
#' @template elipsis_Template
print.summary.hpaSelection <- function (x, ...) 
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  return(print_summary_hpaSelection(x))
}
###
#' Plot hpaSelection random errors approximated density
#' @param x Object of class "hpaSelection"
#' @param y this parameter currently ignored
#' @param is_outcome logical; if TRUE then function plots the graph for outcome equation random errors. 
#' Otherwise plot for selection equation random errors will be plotted.
#' @template elipsis_Template
#' @return This function returns the list containing random error's expected value \code{errors_exp}
#' and variance \code{errors_var} estimates for selection (if \code{is_outcome = TRUE}) or outcome
#' (if \code{is_outcome = FALSE}) equation.
plot.hpaSelection <- function (x, y = NULL, ..., is_outcome = TRUE) 
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  if (!is.null(y))
  {
    warnings("Note that y parameter currently ignored")   
  }
  return(plot_hpaSelection(x, is_outcome))
}
###
#' Calculates AIC for "hpaSelection" object
#' @description This function calculates AIC for "hpaSelection" object
#' @param object Object of class "hpaSelection"
#' @template AIC_Template
#' @template elipsis_Template
AIC.hpaSelection <- function (object, ..., k = 2) 
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  return(AIC_hpaSelection(object, k))
}
###
#' Calculates log-likelihood for "hpaSelection" object
#' @description This function calculates log-likelihood for "hpaSelection" object
#' @usage \method{logLik}{hpaSelection}(object, ...)
#' @param object Object of class "hpaSelection"
#' @template elipsis_Template
logLik.hpaSelection <- function (object, ...) 
{
  if (length(list(...)) > 0)
  {
    warnings("Additional arguments passed throught ... are ignored.")   
  }
  return(logLik_hpaSelection(object))
}