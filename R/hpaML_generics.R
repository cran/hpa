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
    warning("Additional arguments passed through ... are ignored.")   
  }
  return(predict_hpaML(object, as.matrix(newdata)))
}
###
#' Summarizing hpaML Fits
#' @param object Object of class "hpaML"
#' @template elipsis_Template
#' @return This function returns the same list as \code{\link[hpa]{hpaML}} 
#' function changing its class to "summary.hpaML".
summary.hpaML <- function (object, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
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
    warning("Additional arguments passed through ... are ignored.")   
  }
  return(print_summary_hpaML(x))
}

###
#' Calculates log-likelihood for "hpaML" object
#' @description This function calculates log-likelihood for "hpaML" object
#' @usage \method{logLik}{hpaML}(object, ...)
#' @param object Object of class "hpaML"
#' @template elipsis_Template
logLik.hpaML <- function (object, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  lnL <- logLik_hpaML(object)
  attr(lnL, "class") <- "logLik"
  attr(lnL, "df") <- length(as.vector(object$x1))
  
  return(lnL)
}
###
#' Print method for "hpaML" object
#' @param x Object of class "hpaML"
#' @template elipsis_Template
print.hpaML <- function (x, ...)
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
###
#' Plot approximated marginal density using hpaML output
#' @param x Object of class "hpaML"
#' @param y this parameter currently ignored
#' @param ind index of random variable for which
#' approximation to marginal density should be plotted
#' @param given numeric vector of the same length as given_ind
#' from \code{x}. Determines conditional values for the corresponding
#' components. \code{NA} values in \code{given} vector indicate that
#' corresponding random variable is not conditioned. By default all
#' \code{given} components are \code{NA} so unconditional marginal
#' density will be plotted for the \code{ind}-th random variable.
#' @template elipsis_plot_Template
plot.hpaML <- function (x, y = NULL, ..., 
                        ind = 1, given = NULL) 
{
  if (!is.null(y))
  {
    warning("Note that y parameter currently ignored.")   
  }
  
  # Get the data
  my_data <- as.matrix(x$data)
  data_col = ncol(my_data)
  data_row = nrow(my_data)

  # Validate the input
  if (is.null(data_col))
  {
    data_col = 1
  }

  if ((ind < 0) | (ind > data_col))
  {
    stop("data_col argument should be positive integer not greater then number of columns in data.")   
  }
  
  if(is.null(given))
  {
    given <- rep(NA, data_col)
  } else {
    if (length(given) != data_col)
    {
      stop("given argument length should be the same as then number of columns in data.")   
    }
  }
  
  # Get data column corresponding to ind
  h <- NULL
  
  if (data_col == 1)
  {
    h <- matrix(as.numeric(my_data), ncol = 1)
  } else {
    h <- matrix(as.numeric(my_data[, ind]), ncol = 1)
  }
  
  # Assign conditional and marginal indexes
  given_ind_new <- !is.na(given)
  
  omit_ind_new <- rep(TRUE, data_col)
  omit_ind_new[ind] <- FALSE
  omit_ind_new[given_ind_new] <- FALSE
  
  d_ind <- (!given_ind_new) & (!omit_ind_new)
  
  # Estimate moments
  if(is.null(given))
  {
    guven <- vector(mode = "numeric")
  }
  
  h_mean <- ehpa(x = given,
                 pol_coefficients = x$pol_coefficients, 
                 pol_degrees = x$pol_degrees,
                 given_ind = given_ind_new,
                 omit_ind = omit_ind_new,
                 mean = x$mean, sd = x$sd,
                 is_validation = FALSE, 
                 expectation_powers = d_ind)
  h_mean2 <- ehpa(x = given,
                  pol_coefficients = x$pol_coefficients, 
                  pol_degrees = x$pol_degrees,
                  given_ind = given_ind_new,
                  omit_ind = omit_ind_new,
                  mean = x$mean, sd = x$sd,
                  is_validation = FALSE, 
                  expectation_powers = d_ind * 2)
  h_sd <- sqrt(h_mean2 - h_mean ^ 2)
  print(c(h_mean, h_mean2, h_sd, d_ind))
  
  

  # Set maximum and minimum values for plots
  plot_min <- h_mean - 3.8 * h_sd;
  plot_max <- h_mean + 3.8 * h_sd;

  if (!is.na(x$tr_left[ind]))
  {
    plot_min <- max(c(x$tr_left[ind], min(h[, ind])))
  }
  
  if(!is.na(x$tr_right[ind]))
  {
    plot_max <- min(c(x$tr_right[ind], max(h[, ind])))
  }

  # Generate values for plot
  n <- 10000;
  precise <- (plot_max - plot_min) / n;

  x_vec <- matrix(seq(from = plot_min, 
                      by = precise, 
                      to = plot_max), 
                  ncol = 1)
  
  for (i in 1:data_col)
  {
    if (i < ind)
    {
      x_vec <- cbind(given[i], x_vec)
    } else {
      if (i > ind)
      {
        x_vec <- cbind(x_vec, given[i])
      }
    }
  }

  den <- NULL
  
  if(class(x_vec)[1] != "matrix")
  {
    x_vec <- matrix(x_vec, ncol = 1)
  }

  # Calculate density values
  if (any(is.na(x$tr_left)) | any(is.na(x$tr_right)))
  {
    den <- dhpa(x = x_vec,
                pol_coefficients = x$pol_coefficients, 
                pol_degrees = x$pol_degrees,
                given_ind = given_ind_new,
                omit_ind = omit_ind_new,
                mean = x$mean, sd = x$sd,
                is_validation = FALSE)
  } else {
    den <- dtrhpa(x = x_vec,
                  pol_coefficients = x$pol_coefficients, 
                  pol_degrees = x$pol_degrees,
                  given_ind = given_ind_new,
                  omit_ind = omit_ind_new,
                  mean = x$mean, sd = x$sd,
                  tr_left = x$tr_left,
                  tr_right = x$tr_right,
                  is_validation = TRUE)
  }

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
  plot_args$x = x_vec[, ind]
  plot_args$y = den
  
    # make the plot
  do.call(plot, plot_args)
}

#' Extract coefficients from hpaML object
#' @param object Object of class "hpaML"
#' @template elipsis_Template
coef.hpaML <- function (object, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  return(object$x1)
}

#' Extract covariance matrix from hpaML object
#' @param object Object of class "hpaML"
#' @template elipsis_Template
vcov.hpaML <- function (object, ...)
{
  if (length(list(...)) > 0)
  {
    warning("Additional arguments passed through ... are ignored.")   
  }
  
  return(object$cov_mat)
}