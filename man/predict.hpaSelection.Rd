% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hpaSelection_generics.R
\name{predict.hpaSelection}
\alias{predict.hpaSelection}
\title{Predict outcome and selection equation values from hpaSelection model}
\usage{
\method{predict}{hpaSelection}(
  object,
  ...,
  newdata = NULL,
  method = "HPA",
  is_cond = TRUE,
  type = "outcome"
)
}
\arguments{
\item{object}{Object of class "hpaSelection"}

\item{...}{further arguments (currently ignored)}

\item{newdata}{An optional data frame (for \link[hpa]{hpaBinary} and
\link[hpa]{hpaSelection}) or numeric matrix (for \link[hpa]{hpaML})
in which to look for variables with which to predict. If omitted,
the original data frame (matrix) used.}

\item{method}{string value indicating prediction method based 
on hermite polynomial approximation "HPA" or Newey method "Newey".}

\item{is_cond}{logical; if \code{TRUE} (default) then conditional 
predictions will be estimated. Otherwise unconditional 
predictions will be returned.}

\item{type}{character; if "outcome" (default) then predictions for 
selection equation will be estimated according to \code{method}.
If "selection" then selection equation predictions (probabilities) 
will be returned.}
}
\value{
This function returns the list which structure 
depends on \code{method}, \code{is_probit} and \code{is_outcome} values.
}
\description{
This function predicts outcome and selection equation 
values from hpaSelection model.
}
\details{
Note that Newey method can't predict conditional outcomes 
for zero selection equation value. Conditional probabilities for 
selection equation could be estimated only when dependent variable 
from outcome equation is observable.
}
