% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hpaML_generics.R
\name{plot.hpaML}
\alias{plot.hpaML}
\title{Plot approximated marginal density using hpaML output}
\usage{
\method{plot}{hpaML}(x, y = NULL, ..., ind = 1, given = NULL)
}
\arguments{
\item{x}{Object of class "hpaML"}

\item{y}{this parameter currently ignored}

\item{...}{further arguments to be passed to \code{\link[base]{plot}}
function.}

\item{ind}{index of random variable for which
approximation to marginal density should be plotted}

\item{given}{numeric vector of the same length as given_ind
from \code{x}. Determines conditional values for the corresponding
components. \code{NA} values in \code{given} vector indicate that
corresponding random variable is not conditioned. By default all
\code{given} components are \code{NA} so unconditional marginal
density will be plotted for the \code{ind}-th random variable.}
}
\description{
Plot approximated marginal density using hpaML output
}
