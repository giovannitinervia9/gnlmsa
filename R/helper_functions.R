


#-------------------------------------------------------------------------------



#' Numerical Jacobian with respect to parameters
#'
#' Creates a function that computes the numerical Jacobian of a given function
#' with respect to its parameter vector. The returned function evaluates the Jacobian
#' at specified data and parameter values.
#'
#' @param f A function with signature `f(X, theta)` that returns a numeric vector.
#'   X typically represents the data (can be a vector or a matrix),
#'   and theta is the parameter vector for which the Jacobian will be computed.
#'
#' @return A function with signature `function(X, theta)` that computes the Jacobian matrix
#'   of f with respect to theta at the given values. Each row of the resulting Jacobian
#'   corresponds to an element of the output vector, and each column corresponds to
#'   a parameter in theta.
#'
#' @details
#' The function uses the `numDeriv::jacobian()` function to perform numerical differentiation.
#' The returned Jacobian matrix has dimensions `(length(f(X, theta)), length(theta))`.
#'
#' @export
#'
#' @examples
#' # Example 1: Simple exponential model
#' f <- function(X, theta) drop(exp(X %*% theta))
#'
#' # Create a Jacobian function
#' jacobian_fn <- make_jacobian(f)
#'
#' # Compute the Jacobian for iris data
#' X <- as.matrix(iris[, 1:4])
#' theta <- c(0.1, 0.3, -0.2, 0.2)
#' J <- jacobian_fn(X, theta)
#' dim(J)  # Should be (nrow(X), length(theta))
make_jacobian <- function(f) function(X, theta) {
  numDeriv::jacobian(function(theta) f(X, theta), theta)
}


#-------------------------------------------------------------------------------



#' Numerical Hessian with respect to parameters
#'
#' Creates a function that computes the numerical Hessian matrices of a given function
#' with respect to its parameter vector. The returned function evaluates the Hessian
#' at specified data points and parameter values.
#'
#' @param f A function with signature `f(X, theta)` that returns a numeric value.
#'   X typically represents the data (can be a vector or a matrix),
#'   and theta is the parameter vector for which the Hessian will be computed.
#'
#' @return A function with signature `function(X, theta)` that returns a list of Hessian matrices.
#'   Each matrix in the list corresponds to the Hessian of f with respect to theta,
#'   evaluated at a specific data point (row of X).
#'
#' @details
#' The function uses the `numDeriv::hessian` function to perform numerical differentiation.
#' For each observation (row) in X, a separate Hessian matrix is computed, and these
#' matrices are returned as a list. Each Hessian matrix has dimensions `(length(theta), length(theta))`
#' and contains the second-order partial derivatives of f with respect to the parameters.
#'
#' If X has only one row or is a vector, the result will be a list with a single Hessian matrix.
#' If X has multiple rows, the result will be a list of Hessian matrices, one for each row of X.
#'
#' Note that this implementation treats each row of X as a separate data point. The function
#' transposes each row to ensure consistent behavior regardless of whether X is passed as
#' a single vector or as rows in a matrix.
#'
#' @export
#'
#' @examples
#' # Example 1: Simple exponential model
#' f <- function(X, theta) drop(exp(X %*% theta))
#'
#' # Create a Hessian function
#' hessian_fn <- make_hessian(f)
#'
#' # Compute the Hessian for a subset of iris data
#' X <- as.matrix(iris[1:3, 1:4])  # Using only 3 observations for clarity
#' theta <- c(0.1, 0.3, -0.2, 0.2)
#' H_list <- hessian_fn(X, theta)
#'
#' # Examine the result
#' length(H_list)  # Should be 3 (one Hessian per row of X)
#' dim(H_list[[1]])  # Should be (4, 4) (4 parameters)
make_hessian <- function(f) function(X, theta) {
  hi <- function(x) numDeriv::hessian(function(theta) f(t(x), theta), theta)
  apply(as.matrix(X), 1, hi, simplify = FALSE)
}



#-------------------------------------------------------------------------------



#' Maps values between open intervals and the real line
#'
#' This function implements a bijective transformation between an open interval `(lower, upper)`
#' and the real line using a logit-type transformation.
#' @param x a numerical vector of values to map. When `inverse = FALSE`, values must be strictly
#'   between `lower` and `upper`. When `inverse = TRUE`, values can be any real number.
#' @param lower a numerical vector containing the lower bounds of the open intervals.
#' @param upper a numerical vector containing the upper bounds of the open intervals.
#'   Must be greater than the corresponding `lower` values.
#' @param inverse logical indicating whether to map from the open intervals to the real line
#'   (default, `inverse = FALSE`) or from the real line to the open intervals (`inverse = TRUE`).
#' @param derivative logical indicating whether to return the derivative of the mapping function
#'   with respect to x. Only applicable when `inverse = FALSE`.
#'
#' @details
#' Vectors `lower` and `upper` are recycled to match the length of `x` if necessary.
#' The forward mapping (default, `inverse = FALSE`) is \deqn{\log\left(\dfrac{x-\texttt{lower}}{\texttt{upper}-x}\right)}
#' and the inverse mapping (`inverse = TRUE`) is \deqn{\dfrac{\texttt{lower}+\texttt{upper}\exp(x)}{1+\exp(x)}}
#'
#' When `derivative = TRUE`, the function returns the derivative of the forward mapping which is
#' \deqn{\dfrac{\texttt{upper} - \texttt{lower}}{(x - \texttt{lower})(\texttt{upper} - x)}}
#'
#' @return A numerical vector containing the mapped values or the derivative values.
#' @export
#'
#' @examples
#' # Basic usage: map from (0, 1) to R
#' x <- seq(from = 0.1, to = 0.9, length = 10)
#' lower <- 0
#' upper <- 1
#' y <- map_interval(x, lower, upper)
#' y
#'
#' # Map back from R to (0, 1)
#' x_recovered <- map_interval(y, lower, upper, inverse = TRUE)
#' all.equal(x, x_recovered)  # Should be TRUE
#'
#' # Calculate the derivative
#' map_interval(x, lower, upper, derivative = TRUE)
#'
#' # Using different intervals
#' map_interval(c(1.5, 2.5), lower = 1, upper = 3)
#'
#' # Multiple intervals with vectorized lower and upper bounds
#' map_interval(c(0.2, 2.5), lower = c(0, 1), upper = c(1, 3))
map_interval <- function(x, lower, upper, inverse = FALSE, derivative = FALSE) {
  if (any(lower >= upper)) stop("lower must be less than upper")
  if (derivative) {
    (upper - lower)/((x - lower)*(upper - x))
  }
  else if (inverse) {
    (lower + upper*exp(x))/(1 + exp(x))
  } else {
    log((x - lower)/(upper - x))
  }
}



#-------------------------------------------------------------------------------



#' Maps values between half-open intervals and the real line
#'
#' This function implements a bijective transformation between a half-open interval `(lower, Inf)`
#' and the real line using a logarithmic transformation.
#'
#' @param x a numerical vector of values to map. When `inverse = FALSE`, values must be strictly
#'   greater than `lower`. When `inverse = TRUE`, values can be any real number.
#' @param lower a numerical value representing the lower bound of the half-open interval.
#'   Can be any real number (positive, zero, or negative). Default is 0.
#' @param inverse logical indicating whether to map from the half-open interval to the real line
#'   (default, `inverse = FALSE`) or from the real line to the half-open interval (`inverse = TRUE`).
#' @param derivative logical indicating whether to return the derivative of the mapping function
#'   with respect to x. Only applicable when `inverse = FALSE`.
#'
#' @details
#' Despite the function name `map_positive()`, this function can map any half-open interval of the form
#' `(lower, Inf)` to the real line, where `lower` can be any real number.
#'
#' The forward mapping (default, `inverse = FALSE`) is \deqn{\log(x - \texttt{lower})}
#' and the inverse mapping (`inverse = TRUE`) is \deqn{\texttt{lower} + \exp(x)}
#'
#' When `derivative = TRUE`, the function returns the derivative of the forward mapping which is
#' \deqn{\dfrac{1}{x - \texttt{lower}}}
#'
#' @return A numerical vector containing the mapped values or the derivative values.
#' @export
#'
#' @examples
#' # Map from (0, Inf) to R
#' x <- 1:5
#' y <- map_positive(x)
#' y
#'
#' # Map back from R to (0, Inf)
#' x_recovered <- map_positive(y, inverse = TRUE)
#' all.equal(x, x_recovered)  # Should be TRUE
#'
#' # Map from (1, Inf) to R
#' x <- 2:6
#' lower <- 1
#' y <- map_positive(x, lower = lower)
#' y
#'
#' # Map back from R to (1, Inf)
#' x_recovered <- map_positive(y, inverse = TRUE, lower = lower)
#' all.equal(x, x_recovered)  # Should be TRUE
#'
#' # Using a negative lower bound: map from (-5, Inf) to R
#' x <- seq(-4, 10, by = 2)
#' lower <- -5
#' y <- map_positive(x, lower = lower)
#' y
#'
#' # Map back from R to (-5, Inf)
#' x_recovered <- map_positive(y, inverse = TRUE, lower = lower)
#' all.equal(x, x_recovered)  # Should be TRUE
#'
#' # Calculate the derivative for the forward mapping
#' x <- 1:5
#' map_positive(x, derivative = TRUE)
map_positive <- function(x, lower = 0, inverse = FALSE, derivative = FALSE){
  if (derivative) {
    1/(x - lower)
  }
  else if (inverse) {
    lower + exp(x)
  }  else {
    log(x - lower)
  }
}


#-------------------------------------------------------------------------------



#' Maps values between half-open intervals and the real line
#'
#' This function implements a bijective transformation between a half-open interval `(-Inf, upper)`
#' and the real line using a logarithmic transformation.
#'
#' @param x a numerical vector of values to map. When `inverse = FALSE`, values must be strictly
#'   less than `upper`. When `inverse = TRUE`, values can be any real number.
#' @param upper a numerical value representing the upper bound of the half-open interval.
#'   Can be any real number (positive, zero, or negative). Default is 0.
#' @param inverse logical indicating whether to map from the half-open interval to the real line
#'   (default, `inverse = FALSE`) or from the real line to the half-open interval (`inverse = TRUE`).
#' @param derivative logical indicating whether to return the derivative of the mapping function
#'   with respect to x. Only applicable when `inverse = FALSE`.
#'
#' @details
#' This function maps any half-open interval of the form `(-Inf, upper)` to the real line, where
#' `upper` can be any real number.
#'
#' The forward mapping (default, `inverse = FALSE`) is \deqn{\log(\texttt{upper} - x)}
#' and the inverse mapping (`inverse = TRUE`) is \deqn{\texttt{upper} - \exp(x)}
#'
#' When `derivative = TRUE`, the function returns the derivative of the forward mapping which is
#' \deqn{-\dfrac{1}{\texttt{upper} - x}}
#'
#' @return A numerical vector containing the mapped values or the derivative values.
#' @export
#'
#' @examples
#' # Map from (-Inf, 0) to R
#' x <- -(5:1)
#' y <- map_negative(x)
#' y
#'
#' # Map back from R to (-Inf, 0)
#' x_recovered <- map_negative(y, inverse = TRUE)
#' all.equal(x, x_recovered)  # Should be TRUE
#'
#' # Map from (-Inf, 1) to R
#' x <- -(6:2)
#' upper <- 1
#' y <- map_negative(x, upper = upper)
#' y
#'
#' # Map back from R to (-Inf, 1)
#' x_recovered <- map_negative(y, inverse = TRUE, upper = upper)
#' all.equal(x, x_recovered)  # Should be TRUE
#'
#' # Using a negative upper bound: map from (-Inf, -3) to R
#' x <- seq(-10, -4, by = 1)
#' upper <- -3
#' y <- map_negative(x, upper = upper)
#' y
#'
#' # Map back from R to (-Inf, -3)
#' x_recovered <- map_negative(y, inverse = TRUE, upper = upper)
#' all.equal(x, x_recovered)  # Should be TRUE
#'
#' # Calculate the derivative for the forward mapping
#' x <- -(5:1)
#' map_negative(x, derivative = TRUE)
map_negative <- function(x, upper = 0, inverse = FALSE, derivative = FALSE){
  if (derivative) {
    - 1/(upper - x)
  } else if (inverse) {
    upper - exp(x)
  } else {
    log(upper - x)
  }
}



#-------------------------------------------------------------------------------



#' Create mapping functions between constrained parameter spaces and the real line
#'
#' This function generates a set of transformation functions that map parameters from various constrained spaces
#' to the real line and vice versa. It automatically identifies the appropriate transformation for each parameter
#' based on the provided bounds.
#'
#' @param lower a numerical vector containing the lower bounds for each parameter.
#'   Use `-Inf` for unbounded below parameters.
#' @param upper a numerical vector containing the upper bounds for each parameter.
#'   Use `Inf` for unbounded above parameters.
#'
#' @details
#' The function automatically selects the appropriate transformation for each parameter based on its bounds:
#'
#' 1. For parameters with finite bounds `(lower[i], upper[i])`, [map_interval()] is used.
#' 2. For parameters bounded only below `(lower[i], Inf)`, [map_positive()] is used.
#' 3. For parameters bounded only above `(-Inf, upper[i])`, [map_negative()] is used.
#' 4. For completely unbounded parameters `(-Inf, Inf)`, the `identity` function is used.
#'
#' The `map()` function transforms the constrained parameters to unconstrained values on the real line.
#' The `invert()` function performs the inverse transformation from the real line back to the constrained space.
#' The `jacobian()` function returns the derivatives of the transformation for each parameter.
#'
#' @return A list containing three functions:
#'   \describe{
#'     \item{`map(par)`}{A function that maps parameters from the constrained space to the real line.}
#'     \item{`invert(par)`}{A function that maps parameters from the real line back to the constrained space.}
#'     \item{`jacobian(par)`}{A function that returns the derivative of the transformation for each parameter.}
#'   }
#'   All three functions expect a numeric vector `par` of the same length as `lower` and `upper`.
#'
#' @export
#'
#' @examples
#' # Create mapping functions for parameters with different constraints
#' lower <- c(-Inf, 0, -1)
#' upper <- c(Inf, Inf, 1)
#'
#' # Interpretation:
#' # - First parameter: unbounded (-Inf, Inf)
#' # - Second parameter: bounded below (0, Inf)
#' # - Third parameter: bounded on both sides (-1, 1)
#'
#' map_functions <- make_map_function(lower, upper)
#'
#' # Extract the mapping functions
#' map <- map_functions$map
#' invert <- map_functions$invert
#' jacobian <- map_functions$jacobian
#'
#' # Define parameter values in the constrained space
#' x <- c(0, 3, 0.2)
#'
#' # Map to unconstrained space
#' y <- map(x)
#' y
#'
#' # Map back to constrained space
#' x_recovered <- invert(y)
#' x_recovered
#'
#' # Verify roundtrip accuracy
#' x - x_recovered
#' all.equal(x, x_recovered)
#'
#' # Calculate the Jacobian at point x
#' j <- jacobian(x)
#' j
make_map_function <- function(lower, upper){
  if (length(lower) != length(upper)) {
    stop("lower and upper must be of the same length")
  }

  if (any(lower >= upper)) stop("lower must be less than upper")

  npar <- length(lower)
  f_string <- character(npar)
  f_string_inverse <- f_string
  j_string <- f_string
  for(i in 1:npar){
    if (is.infinite(lower[i]) & is.infinite(upper[i])) {
      f_string[i] <- f_string_inverse[i] <- paste0("identity(par[", deparse(i), "])")
      j_string[i] <- "1"
    }
    else if (!is.infinite(lower[i]) & is.infinite(upper[i])) {
      f_string[i] <- paste0("map_positive(par[", deparse(i), "], inverse = FALSE, lower = ", deparse(lower[i]), ")")
      f_string_inverse[i] <- paste0("map_positive(par[", deparse(i), "], inverse = TRUE, lower = ", deparse(lower[i]), ")")
      j_string[i] <- paste0("map_positive(par[", deparse(i), "], inverse = FALSE, derivative = TRUE, lower = ", deparse(lower[i]), ")")
    }
    else if (is.infinite(lower[i]) & !is.infinite(upper[i])) {
      f_string[i] <- paste0("map_negative(par[", deparse(i), "], inverse = FALSE, upper = ", deparse(upper[i]), ")")
      f_string_inverse[i] <- paste0("map_negative(par[", deparse(i), "], inverse = TRUE, upper = ", deparse(upper[i]), ")")
      j_string[i] <- paste0("map_negative(par[", deparse(i), "], inverse = FALSE, derivative = TRUE, upper = ", deparse(upper[i]), ")")
    }
    else if(!is.infinite(lower[i]) & !is.infinite(upper[i])){
      f_string[i] <- paste0("map_interval(par[", deparse(i), "], lower = ", deparse(lower[i]), ", upper = ", deparse(upper[i]), ", inverse = FALSE)")
      f_string_inverse[i] <- paste0("map_interval(par[", deparse(i), "], lower = ", deparse(lower[i]), ", upper = ", deparse(upper[i]), ", inverse = TRUE)")
      j_string[i] <- paste0("map_interval(par[", deparse(i), "], lower = ", deparse(lower[i]), ", upper = ", deparse(upper[i]), ", inverse = FALSE, derivative = TRUE)")
    }
  }
  map <- function(par) {
    code <- paste0("c(", paste(f_string, collapse = ", "), ")")
    eval(parse(text = code))
  }
  invert <- function(par) {
    code <- paste0("c(", paste(f_string_inverse, collapse = ", "), ")")
    eval(parse(text = code))
  }
  jacobian <- function(par) {
    code <- paste0("c(", paste(j_string, collapse = ", "), ")")
    eval(parse(text = code))
  }
  list(map = map, invert = invert, jacobian = jacobian)
}



#-------------------------------------------------------------------------------



#' Sample new parameter values
#'
#' This function samples new parameter values using a multivariate normal distribution,
#' based on previous values and a variance-covariance matrix. If sampling with mvtnorm::rmvnorm fails,
#' it falls back to independent sampling using stats::rnorm.
#'
#' @param par Numeric vector containing the previous parameter values.
#' @param v Square symmetric variance-covariance matrix used to sample errors.
#' @param mult Numeric vector to multiply the sampled errors. Can be a single value
#'   (applied to all parameters) or a vector of the same length as `par`. Default is 1.
#' @param npar Number of parameters. Default is `length(par)`.
#'
#' @return Numeric vector containing the new parameter values.
#' @export
#' @importFrom stats rnorm
#' @importFrom mvtnorm rmvnorm
#'
#' @examples
#' # Define initial parameters and their bounds
#' par <- c(0, 2, .2)
#' lower <- c(-Inf, 0, -1)
#' upper <- c(Inf, Inf, 1)
#'
#' # Create mapping functions to handle parameter bounds
#' map_functions <- make_map_function(lower, upper)
#' map <- map_functions$map
#' invert <- map_functions$invert
#' jacobian <- map_functions$jacobian
#'
#' # Map parameters to unbounded space
#' mapped_pars <- map(par)
#'
#' # Calculate Jacobian adjustment for covariance matrix based on delta method
#' j <- diag(jacobian(par))
#' v <- diag(1e-06, length(par))
#' mapped_v <- j %*% v %*% t(j)
#'
#' # Define step size multipliers for different parameters
#' mult <- c(3, 1, .2)
#'
#' # Sample new parameters in unbounded space
#' new_par_mapped <- sample_par(par = mapped_pars, v = mapped_v, mult = mult)
#'
#' # Map parameters back to original constrained space
#' invert(new_par_mapped)
#'
sample_par <- function(par, v, mult = 1, npar = length(par)){
  if(nrow(v) != ncol(v)) stop("v must be a square matrix")
  if(!isSymmetric(v)) stop("v must be a symmetric matrix")
  if(nrow(v) != npar) stop("dimension of variance-covariance matrix and par must coincide")
  if(length(mult) != 1 & length(mult) != npar) stop("mult must be of dimension 1 or npar")
  eps <- tryCatch(c(mvtnorm::rmvnorm(1, sigma = v)),
                  error = function(e) {
                    stats::rnorm(npar, mean = 0, sd = sqrt(diag(v)))
                  }
  )
  par + mult*eps
}














