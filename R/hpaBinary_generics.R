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
      warning("Additional arguments passed through ... are ignored.")   
  }
  return(predict_hpaBinary(object, newdata, is_prob))
}
###
#' Summarizing hpaBinary Fits
#' @param object Object of class "hpaBinary"
#' @template elipsis_Template
#' @return This function returns the same list as \code{\link[hpa]{hpaBinary}} 
#' function changing its class to "summary.hpaBinary".
summary.hpaBinary <- function (object, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  return(summary_hpaBinary(object))
}
###
#' Summary for "hpaBinary" object
#' @param x Object of class "hpaBinary"
#' @template elipsis_Template
print.summary.hpaBinary <- function (x, ...) 
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  return(print_summary_hpaBinary(x))
}
###
#' Plot hpaBinary random errors approximated density
#' @param x Object of class "hpaBinary"
#' @param y this parameter currently ignored
#' @template elipsis_plot_Template
plot.hpaBinary <- function (x, y = NULL, ...) 
{
  if (!is.null(y))
  {
    warning("Note that y parameter currently ignored.")   
  }
  
  # Load data from the model
  pol_coefficients <- as.numeric(unlist(x["pol_coefficients"]))
  pol_degrees <- as.numeric(unlist(x["pol_degrees"]))

  mean <- as.numeric(x["mean"])
  sd <- as.numeric(x["sd"])
  
  # Adjust precision
  errors_exp <- as.numeric(x["errors_exp"])
  errors_var <- as.numeric(x["errors_var"])

  plot_min <- errors_exp - 3.8 * sqrt(errors_var);
  plot_max <- errors_exp + 3.8 * sqrt(errors_var);
  
  n <- 10000;
  
  precise <- (plot_max - plot_min) / n;
  
  x_vec <- plot_min + cumsum(rep(precise, n))

  den <- dhpa(x = x_vec,
              pol_coefficients = pol_coefficients, 
              pol_degrees = pol_degrees,
              mean = mean, sd = sd)
  
  den_min = min(den)
  den_max = max(den)
  
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
#' Calculates log-likelihood for "hpaBinary" object
#' @description This function calculates log-likelihood for "hpaBinary" object
#' @usage \method{logLik}{hpaBinary}(object, ...)
#' @param object Object of class "hpaBinary"
#' @template elipsis_Template
logLik.hpaBinary <- function (object, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  lnL <- logLik_hpaBinary(object)
  attr(lnL, "class") <- "logLik"
  attr(lnL, "df") <- length(as.vector(object$x1))
  
  return(lnL)
}
###
#' Print method for "hpaBinary" object
#' @param x Object of class "hpaBinary"
#' @template elipsis_Template
print.hpaBinary <- function (x, ...) 
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
  cat(paste("Log-likelihood function value is:", round(x$'log-likelihood', 3), "\n"))
  cat("---\n")
  cat("Please, use summary() function to get additional information\n")
}

#' Extract coefficients from hpaBinary object
#' @param object Object of class "hpaBinary"
#' @template elipsis_Template
coef.hpaBinary <- function (object, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  return(object$x1)
}

#' Extract covariance matrix from hpaBinary object
#' @param object Object of class "hpaBinary"
#' @template elipsis_Template
vcov.hpaBinary <- function (object, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  return(object$cov_mat)
}