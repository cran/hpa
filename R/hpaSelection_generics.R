#' Predict outcome and selection equation values from hpaSelection model
#' @description This function predicts outcome and selection equation 
#' values from hpaSelection model.
#' @param object Object of class "hpaSelection"
#' @param method string value indicating prediction method based 
#' on hermite polynomial approximation "HPA" or Newey method "Newey".
#' @template newdata_Template
#' @param is_cond logical; if \code{TRUE} (default) then conditional 
#' predictions will be estimated. Otherwise unconditional 
#' predictions will be returned.
#' @param type character; if "outcome" (default) then predictions for 
#' selection equation will be estimated according to \code{method}.
#' If "selection" then selection equation predictions (probabilities) 
#' will be returned.
#' @template elipsis_Template
#' @details Note that Newey method can't predict conditional outcomes 
#' for zero selection equation value. Conditional probabilities for 
#' selection equation could be estimated only when dependent variable 
#' from outcome equation is observable.
#' @return This function returns the list which structure 
#' depends on \code{method}, \code{is_probit} and \code{is_outcome} values.
predict.hpaSelection <- function (object,  ...,
                                  newdata = NULL, 
                                  method = "HPA", 
                                  is_cond = TRUE, 
                                  type = "outcome")
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  is_outcome <- TRUE
  if(type == "selection")
  {
    is_outcome <- FALSE
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
#' @return This function returns the same list 
#' as \code{\link[hpa]{hpaSelection}} 
#' function changing its class to "summary.hpaSelection".
summary.hpaSelection <- function (object, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  return(summary_hpaSelection(object))
}
###
#' Summary for "hpaSelection" object
#' @param x Object of class "hpaSelection"
#' @template elipsis_Template
print.summary.hpaSelection <- function (x, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  return(print_summary_hpaSelection(x))
}
###
#' Plot hpaSelection random errors approximated density
#' @param x Object of class "hpaSelection"
#' @param y this parameter currently ignored
#' @param type character; if "outcome" then function plots the graph for 
#' outcome equation random errors, if "selection" then plot for selection 
#' equation random errors will be generated.
#' @template elipsis_plot_Template
plot.hpaSelection <- function (x, y = NULL, ..., type = "outcome") 
{
  if (!is.null(y))
  {
    warning("Note that y parameter currently ignored.")   
  }
  
  if(!(type %in% c("outcome", "selection")))
  {
    stop("type argument should be either outcome or selection")
  }
  
  is_outcome <- TRUE
  if(type == "selection")
  {
    is_outcome <- FALSE
  }
  
  # load data from the model
  
  selection_mean <- as.numeric(x["selection_mean"])
  outcome_mean <- as.numeric(x["outcome_mean"])
  selection_sd <- as.numeric(x["selection_sd"])
  outcome_sd <- as.numeric(x["outcome_sd"])
  
  mean = c(selection_mean, outcome_mean)
  sd = c(selection_sd, outcome_sd)
  
  pol_degrees <- as.numeric(unlist(x["pol_degrees"]))
  pol_coefficients <- as.numeric(unlist(x["pol_coefficients"]))
  
  errors_exp <- ehpa(pol_coefficients = pol_coefficients,
                     pol_degrees = pol_degrees,
                     mean = mean, sd = sd,
                     expectation_powers = c(1 - is_outcome, is_outcome),
                     is_validation = FALSE)
  
  errors_exp_2 <- ehpa(pol_coefficients = pol_coefficients,
                       pol_degrees = pol_degrees,
                       mean = mean, sd = sd,
                       expectation_powers = 2 * c((1 - is_outcome), is_outcome),
                       is_validation =  FALSE)
  
  errors_var <- errors_exp_2 - (errors_exp * errors_exp)

  # Adjust precision
  
  plot_min <- errors_exp - 3.8 * sqrt(errors_var);
  plot_max <- errors_exp + 3.8 * sqrt(errors_var);

  n <- 10000;
  
  precise <- (plot_max - plot_min) / n

  x_vec <- plot_min + cumsum(rep(precise, n))

  # calculate densities
  
  den <- dhpa(cbind(x_vec, x_vec),
              pol_coefficients = pol_coefficients,
              pol_degrees = pol_degrees,
              omit_ind = c(is_outcome, !is_outcome),
              mean = mean, sd = sd, 
              is_validation = FALSE);

  den_min <- min(den);
  den_max <- max(den);
  
  # Plot the result
  
    # prepare the arguments
  plot_args <- list(...)
  plot_args_names <- names(plot_args)
  if(!("xlim" %in% plot_args_names))
  {
    plot_args$xlim <- c(plot_min, plot_max)
  }
  if(!("type" %in% plot_args_names))
  {
    plot_args$type <- "l"
  }
  if(!("lwd" %in% plot_args_names))
  {
    plot_args$lwd <- 3
  }
  if(!("main" %in% plot_args_names))
  {
    plot_args$main <- "Random Errors Density Approximation Plot"
  }
  if(!("xlab" %in% plot_args_names))
  {
    plot_args$xlab <- "value"
  }
  if(!("ylab" %in% plot_args_names))
  {
    plot_args$ylab <- "density"
  }
  plot_args$x = x_vec
  plot_args$y = den
  
    # make the plot
  do.call(plot, plot_args)
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
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  lnL <- logLik_hpaSelection(object)
  attr(lnL, "class") <- "logLik"
  attr(lnL, "df") <- length(as.vector(object$x1))
  
  return(lnL)
}
###
#' Print method for "hpaSelection" object
#' @param x Object of class "hpaSelection"
#' @template elipsis_Template
print.hpaSelection <- function (x, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  cat(paste("It is the object of class",class(x),"\n"))
  cat("It contains the following elements:\n")
  cat(names(x), sep = ", ")
  cat("\n")
  cat("---\n")
  cat("Estimation results:\n")
  print(x$results)
  cat("---\n")
  cat(paste("Log-likelihood function value is:", 
            round(x$'log-likelihood', 3), "\n"))
  cat("---\n")
  cat("Please, use summary() function to get additional information\n")
}

#' Extract coefficients from hpaSelection object
#' @param object Object of class "hpaSelection"
#' @param type character; if "all" (default) then all estimated parameters
#' values will be returned. If "selection" then selection equation coefficients
#' estimates will be provided. If "outcome" then outcome equation coefficients
#' estimates will be returned.
#' @template elipsis_Template
coef.hpaSelection <- function (object, ..., type = "all")
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  if (type == "outcome")
  {
    return(object$outcome_coef)
  }
  
  if (type == "selection")
  {
    return(object$selection_coef)
  }
  
  return (object$x1)
}

#' Extract covariance matrix from hpaSelection object
#' @param object Object of class "hpaSelection"
#' @template elipsis_Template
vcov.hpaSelection <- function (object, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  return(object$cov_mat)
}