% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{reparametrize}
\alias{reparametrize}
\title{Reparameterize Gradients and Hessians for \code{gnlmsa_fit}}
\usage{
reparametrize(theta, g, h, map_functions)
}
\arguments{
\item{theta}{Numeric vector of the parameters in the constrained space.}

\item{g}{Numeric vector representing the gradient of the log-likelihood function in the constrained space.}

\item{h}{Numeric matrix representing the hessian of the log-likelihood function in the constrained space.}

\item{map_functions}{a list as created by \code{\link[=make_map_function]{make_map_function()}}.}
}
\value{
A named list with two elements:
\describe{
\item{\code{g_map}}{Numeric vector representing the gradient of the log-likelihood function in the unconstrained space.}
\item{\code{h_map}}{Numeric matrix representing the hessian of the log-likelihood function in the unconstrained space.}
}
}
\description{
In the \pkg{gnlmsa} framework, model parameters may be
constrained (e.g., positive‐only, bounded between 0 and 1).
\pkg{gnlmsa_fit} performs likelihood optimization via Newton-Raphson algorithm on the
unconstrained scale (see \code{\link[=make_map_function]{make_map_function()}} for details about the transformation).
\code{reparametrize()} automates the chain‐rule adjustment
of both the gradient and Hessian of the log-likelihood
when moving from the constrained space to the unconstrained space.
}
\keyword{internal}
