% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gnlmsa_fit.R
\name{nr_control}
\alias{nr_control}
\title{Control Parameters for Newton-Raphson Algorithm}
\usage{
nr_control(
  maxit = 100,
  tol = 1e-06,
  expected = FALSE,
  unconstrained = TRUE,
  verbose = TRUE,
  regularization = 0.001,
  max_step_halving = 1000,
  max_decr = 10
)
}
\arguments{
\item{maxit}{Integer. Maximum number of iterations.}

\item{tol}{Numeric. Convergence tolerance on relative log-likelihood increase.}

\item{expected}{Logical; if \code{TRUE}, the expected Fisher information is used in the optimization process.
If \code{FALSE} (default), the observed Hessian is used.}

\item{unconstrained}{Logical; if \code{TRUE}, the optimization is performed in the unconstrained space (for detail see \code{\link[=make_map_function]{make_map_function()}} and \code{\link[=reparametrize]{reparametrize()}}.
If \code{FALSE} (default), the observed Hessian is the original parameter space.}

\item{verbose}{Logical; if \code{TRUE} (default), progress and warnings are printed.}

\item{regularization}{Numeric. A (small) constant added to the diagonal of the Hessian if it results non invertible.}

\item{max_step_halving}{Integer. Maximum number of step halving iteration if the log-likelihood decrease within the optimization iteration.}

\item{max_decr}{Integer. Maximum number of iterations with decreasing log-likelihood before stopping the algorithm.}
}
\value{
A list containing the control parameters to be used in the Newton-Raphson routine.
}
\description{
Defines a set of user-configurable options to control the behavior of the Newton-Raphson
routine used in \code{gnlmsa} for optimization of the log-likelihood Generalized Non-Linear Models.
}
