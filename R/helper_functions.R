


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
#' @return A function with signature function(X, theta) that computes the Jacobian matrix
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
#' @param f A function with signature f(X, theta) that returns a numeric value.
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



#' Maps open intervals to real line
#'
#' @param x a numerical vector of values to map to the real line or to take bake to the open intervals.
#' @param lower a numerical vector containing the lower bounds of the open intervals.
#' @param upper a numerical vector containing the upper bounds of the open intervals.
#' @param inverse logical indicating whether to map from the open intervals to the real line (default, `inverse = FALSE`) or whether to map from the real line to the open intervals.
#' @param derivative logical indicating whether to return the derivative of the map function.
#'
#' @return a numerical vector containing the mapped values or the derivative.
#' @export
#'
#' @examples
#' # Create a numerical vector
#' x <- seq(from = 0.01, to = 0.99, length = 10)
#' lower <- 0
#' upper <- 1
#'
#' # map from (0, 1) to R
#' y <- map_interval(x, lower, upper)
#'
#' # map back from R to (0, 1)
#' map_interval(y, lower, upper, inverse = TRUE)
#'
#' # derivative
#' map_interval(x, lower, upper, derivative = TRUE)
map_interval <- function(x, lower, upper, inverse = FALSE, derivative = FALSE){
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



#' Map the positive line to the real line
#'
#' @param x a numerical vector of values to map to the real line or to take back to the positive line.
#' @param inverse logical indicating whether to map from the positive line to the real line (default, `inverse = FALSE`) or whether to map from the real line to the positive line (`inverse = TRUE`).
#' @param lower value of the lower bound in order to map from \eqn{(\texttt{lower}, \infty) \rightarrow \mathbb{R}}. `lower` could also be a negative value. Default is 0.
#' @param derivative logical indicating whether to return the derivative of the map function.
#'
#' @return a numerical vector containing the mapped values or the derivative.
#' @export
#'
#' @examples
#' # (0, inf) to R
#' x <- 1:5
#' y <- map_positive(x)
#' map_positive(y, inverse = TRUE)
#'
#' # (a, Inf) to R
#'
#' x <- 2:6
#' lower <- 1
#' y <- map_positive(x, lower = lower)
#' map_positive(y, inverse = TRUE, lower = lower)
#'
#' # derivative
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



#' Map the negative line to the real line
#'
#' @param x a numerical vector of values to map to the real line or to take back to the negative line.
#' @param inverse logical indicating whether to map from the negative line to the real line (default, `inverse = FALSE`) or whether to map from the real line to the negative line (`inverse = TRUE`).
#' @param upper value of the upper bound in order to map from \eqn{(-\infty, \texttt{upper}) \rightarrow \mathbb{R}}. `upper` could also be a positive value. Default is 0.
#' @param derivative logical indicating whether to return the derivative of the map function.
#'
#' @return a numerical vector containing the mapped values or the derivative.
#' @export
#'
#' @examples
#' # (-Inf, 0) to R
#' x <- -(5:1)
#' y <- map_negative(x)
#' map_negative(y, inverse = TRUE)
#'
#' # (-Inf, a) to R
#' x <- -(6:2)
#' upper <- 1
#' y <- map_negative(x, upper = upper)
#' map_negative(y, inverse = TRUE, upper = upper)
#'
#' # derivative
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



#' Helper function to create a function to map parameters from constrained parameter spaces to the real line
#'
#' @param lower a numerical vector containing the lower bounds. There must be an element for each parameter of the model.
#' @param upper a numerical vector containing the upper bounds. There must be an element for each parameter of the model.
#'
#' @return a list with three functions:
#'  - `map()`: function that maps from the constrained space to the real line.
#'  - `invert()`: inverse of `map()` that maps back from the real line to the constrained space.
#'  - `jacobian()`: jacobian of `map()`.
#'
#' @export
#'
#' @examples
#'
#' lower <- c(-Inf, 0, -1)
#' upper <- c(Inf, Inf, 1)
#'
#' map_functions <- make_map_function(lower, upper)
#'
#' map <- map_functions$map
#' invert <- map_functions$invert
#' jacobian <- map_functions$jacobian
#'
#' x <- c(0, 3, 0.2)
#' map(x)
#' invert(map(x))
#' x - invert(map(x))
#' jacobian(x)
make_map_function <- function(lower, upper){

  if (length(lower) != length(upper)) {
    stop("lower and upper must be of the same length")
  }

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
      f_string[i] <- paste0("map_negative(par[", deparse(i), "], inverse = FALSE, upper = ", deparse(lower[i]), ")")
      f_string_inverse[i] <- paste0("map_negative(par[", deparse(i), "], inverse = TRUE, upper = ", deparse(lower[i]), ")")
      j_string[i] <- paste0("map_negative(par[", deparse(i), "], inverse = FALSE, derivative = TRUE, upper = ", deparse(lower[i]), ")")
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



