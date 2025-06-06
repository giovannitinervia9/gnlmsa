% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{map_positive}
\alias{map_positive}
\title{Create mapping functions between half-open intervals and the real line}
\usage{
map_positive()
}
\value{
A list of six functions: \code{map}, \code{invert}, \code{map_jacobian},
\code{map_hessian}, \code{invert_jacobian}, \code{invert_hessian}.
}
\description{
Returns a list of functions implementing a bijective transformation between
a half-open interval \verb{(lower, Inf)} and the real line using a logarithmic transformation,
together with their derivatives and inverses.
}
\details{
The returned list contains:
\describe{
\item{map(x, lower)}{Maps from (lower, Inf) to the real line: \deqn{f(x) = \log(x - \texttt{lower})}}
\item{invert(x, lower)}{Inverse mapping: \deqn{f^{-1}(z) = \texttt{lower} + \exp(z)}}
\item{map_jacobian(x, lower)}{First derivative of the forward map: \deqn{f'(x) = 1/(x - \texttt{lower})}}
\item{map_hessian(x, lower)}{Second derivative of the forward map: \deqn{f''(x) = -1/(x - \texttt{lower})^2}}
\item{invert_jacobian(x, lower)}{First derivative of the inverse map: \deqn{(f^{-1})'(z) = \exp(z)}}
\item{invert_hessian(x, lower)}{Second derivative of the inverse map: \deqn{(f^{-1})''(z) = \exp(z)}}
}
}
\examples{
f <- map_positive()
x <- 1:5
f$map(x, lower = 0)
f$invert(f$map(x, lower = 0), lower = 0)
f$map_jacobian(x, lower = 0)
f$invert_jacobian(f$map(x, lower = 0), lower = 0)
f$map_hessian(x, lower = 0)
f$invert_hessian(f$map(x, lower = 0), lower = 0)
}
