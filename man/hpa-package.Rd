\name{hpa-package}
\alias{hpa-package}
\alias{hpa}
\docType{package}
\title{
  Provides Tools for Multivariate Distributions Hermite Polynomial Approximation
}
\description{
  This package provides functions to calculate conditional marginal densities, moments, cumulative distribution functions as well as binary choice and sample selection models based on hermite polynomial approximation.}

\details{
Densities hermite polynomial approximation approach has been proposed by A. Gallant and D. W. Nychka in 1987. The main idea is to approximate unknown distribution density with hermite polynomial. In this framework hermite polynomial represents adjusted (to insure integration to 1) product of squared polynomial and normal distribution densities. This package provides functions allowing to use this approach and some related semi-nonparametric models.

Functions starting with single letter and ending on "hpa" provide hermite polynomial approximations itself. The letters are "d", "p", "i" and "e". Letter "d" states for joint density so function \code{\link[hpa]{dhpa}} allows to approximate joint densities with hermite polynomials. Similary functions starting with letters "p", "i", and "e" indicate methods that allow to calculate cumulative distribution function (\code{\link[hpa]{phpa}}), interaval distribution function (\code{\link[hpa]{ihpa}}) and expected value of random's vector components product (\code{\link[hpa]{ehpa}}) approximations. If these letters follow "tr" postfix it means that this functions provide truncated distribution approximation tools. For example, function \code{\link[hpa]{dtrhpa}} allows to approximate truncated joint densities with hermite polynomials. The other two truncated distributions approximating functions are \code{\link[hpa]{itrhpa}} and \code{\link[hpa]{etrhpa}}. Note that all these functions (except \code{\link[hpa]{etrhpa}}) are vectorized respect to some of their arguments and allow to approximate conditional and marginal distributions too.

In order to perform semi-nonparametric maximum likelihood estimation based on hermite polynomial approximation please use function \code{\link[hpa]{hpaML}}. Functions \code{\link[hpa]{hpaBinary}} and \code{\link[hpa]{hpaSelection}} allows to estimate semi-nonparametric single index binary choice and sample selection models correspondlingly.

This package also provides functions allowing to estimate normal distribution  and truncated normal distribution moments. These functions are \code{\link[hpa]{normalMoment}} and \code{\link[hpa]{truncatedNormalMoment}}.

This package has been written completely in \pkg{Rcpp} in order to insure high speed of calculations. However sometimes (where it has no sufficient impact upon calculation speed) calls to the functions from R environments \code{stats}, \code{base} and \code{graphics} have been made.

}
\author{
Potanin Bogdan, bogdanpotanin@gmail.com.

Maintainer: Potanin Bogdan <bogdanpotanin@gmail.com>
}
\references{
A. Gallant and D. W. Nychka (1987) <doi:10.2307/1913241>
}
\keyword{ package }
\seealso{
  \code{\link[hpa]{polynomialIndex}}
}
\examples{
  \donttest{
##Let's approximate some three random variables joint density function at point (0,1, 0.2, 0.3)
##with hermite polynomial of (1,2,3) degrees which polynomial coefficients equals 1 except
##coefficient related to x1*(x^3) polynomial element which equals 2.
##Also suppose that normal density related mean vector equals (1.1, 1.2, 1.3) while
##standard deviations vector is (2.1, 2.2, 2.3).

#Prepare initial values
x <- matrix(c(0.1, 0.2, 0.3), nrow=1)
mean <- c(1.1, 1.2, 1.3)
sd <- c(2.1, 2.2, 2.3)
pol_degrees <- c(1, 2, 3)

#Create polynomial powers and indexes correspondence matrix
pol_ind <- polynomialIndex(pol_degrees)

#Set all polynomial coefficients to 1
pol_coefficients <- rep(1, ncol(pol_ind))
pol_degrees_n <- length(pol_degrees)

#Assign coefficient 2 to the polynomial element(x1 ^ 1)*(x2 ^ 0)*(x3 ^ 2)
pol_coefficients[which(colSums(pol_ind == c(1, 0, 2)) == pol_degrees_n)] <- 2
#Visualize correspondence between polynomial elements and their coefficients
as.data.frame(rbind(pol_ind, pol_coefficients),
	row.names = c("x1 power", "x2 power", "x3 power", "coefficients"),
	optional = TRUE)
printPolynomial(pol_degrees, pol_coefficients)

#Calculate density approximation at point x
dhpa(x = x,
	pol_coefficients = pol_coefficients, pol_degrees = pol_degrees,
	mean = mean, sd = sd)
	
#Condition second component to be 0.5
#Substitute x second component with conditional value 0.5
x <- matrix(c(0.1, 0.5, 0.3), nrow = 1)

#Set TRUE to the second component indicating that it is conditioned
given_ind <- c(FALSE, TRUE, FALSE)

#Calculate conditional (on x2=0.5) density approximation at point x
dhpa(x = x,
	pol_coefficients = pol_coefficients, pol_degrees = pol_degrees,
	mean = mean, sd = sd,
	given_ind = given_ind)
	
#Consider third component marginal distribution
#conditioned on the second component 0.5 value
#Set TRUE to the first component indicating that it is omitted
omit_ind <- c(TRUE, FALSE, FALSE)

#Calculate conditional (on x2=0.5) marginal (for x3) density approximation at point x
dhpa(x = x,
	pol_coefficients = pol_coefficients, pol_degrees = pol_degrees,
	mean = mean, sd = sd,
	given_ind = given_ind, omit_ind = omit_ind) 
  }
}
