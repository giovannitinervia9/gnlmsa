


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



#' Create mapping functions between open intervals and the real line
#'
#' Returns a list of functions implementing a bijective transformation between
#' an open interval `(lower, upper)` and the real line, together with their derivatives and inverses.
#'
#' @details
#' The returned list contains:
#' \describe{
#'   \item{map(x, lower, upper)}{Maps from (lower, upper) to the real line: \deqn{f(x) = \log\left(\dfrac{x - \texttt{lower}}{\texttt{upper} - x}\right)}}
#'   \item{invert(x, lower, upper)}{Inverse mapping: \deqn{f^{-1}(z) = \dfrac{\texttt{lower} + \texttt{upper} \cdot e^z}{1 + e^z}}}
#'   \item{map_jacobian(x, lower, upper)}{First derivative of the forward map: \deqn{f'(x) = \dfrac{\texttt{upper} - \texttt{lower}}{(x - \texttt{lower})(\texttt{upper} - x)}}}
#'   \item{map_hessian(x, lower, upper)}{Second derivative of the forward map: \deqn{f''(x) = -\dfrac{(\texttt{lower} - \texttt{upper})(\texttt{lower} + \texttt{upper} - 2x)}{(x - \texttt{lower})^2 (\texttt{upper} - x)^2}}}
#'   \item{invert_jacobian(x, lower, upper)}{First derivative of the inverse map: \deqn{(f^{-1})'(z) = \dfrac{e^z(\texttt{upper} - \texttt{lower})}{\left(1 + e^z\right)^2}}}
#'   \item{invert_hessian(x, lower, upper)}{Second derivative of the inverse map: \deqn{(f^{-1})''(z) = \dfrac{e^z(e^z - 1)(\texttt{lower} - \texttt{upper})}{\left(1 + e^z\right)^3}}}
#' }
#'
#' @return A list of six functions: \code{map}, \code{invert}, \code{map_jacobian},
#'   \code{map_hessian}, \code{invert_jacobian}, \code{invert_hessian}.
#'
#' @export
#'
#' @examples
#' f <- map_interval()
#' x <- seq(0.1, 0.9, length.out = 10)
#' f$map(x, lower = 0, upper = 1)
#' f$invert(f$map(x, lower = 0, upper = 1), lower = 0, upper = 1)
#' f$map_jacobian(x, lower = 0, upper = 1)
#' f$invert_jacobian(f$map(x, lower = 0, upper = 1), lower = 0, upper = 1)
#' f$map_hessian(x, lower = 0, upper = 1)
#' f$invert_hessian(f$map(x, lower = 0, upper = 1), lower = 0, upper = 1)
map_interval <- function() {

  map <- function(x, lower, upper) {
    log((x - lower) / (upper - x))
  }

  invert <- function(x, lower, upper) {
    (lower + upper * exp(x)) / (1 + exp(x))
  }

  map_jacobian <- function(x, lower, upper) {
    (upper - lower) / ((x - lower) * (upper - x))
  }

  map_hessian <- function(x, lower, upper) {
    ((lower - upper)*(lower + upper - 2*x))/(((x - lower)^2)*((upper - x)^2))
  }

  invert_jacobian <- function(x, lower, upper) {
    (exp(x)*(upper - lower))/(1 + exp(x))^2
  }

  invert_hessian <- function(x, lower, upper) {
    (exp(x)*(exp(x) - 1)*(lower - upper))/(1 + exp(x))^3
  }

  list(map = map, invert = invert,
       map_jacobian = map_jacobian, map_hessian = map_hessian,
       invert_jacobian = invert_jacobian, invert_hessian = invert_hessian)
}



#-------------------------------------------------------------------------------



#' Create mapping functions between half-open intervals and the real line
#'
#' Returns a list of functions implementing a bijective transformation between
#' a half-open interval `(lower, Inf)` and the real line using a logarithmic transformation,
#' together with their derivatives and inverses.
#'
#' @details
#' The returned list contains:
#' \describe{
#'   \item{map(x, lower)}{Maps from (lower, Inf) to the real line: \deqn{f(x) = \log(x - \texttt{lower})}}
#'   \item{invert(x, lower)}{Inverse mapping: \deqn{f^{-1}(z) = \texttt{lower} + \exp(z)}}
#'   \item{map_jacobian(x, lower)}{First derivative of the forward map: \deqn{f'(x) = 1/(x - \texttt{lower})}}
#'   \item{map_hessian(x, lower)}{Second derivative of the forward map: \deqn{f''(x) = -1/(x - \texttt{lower})^2}}
#'   \item{invert_jacobian(x, lower)}{First derivative of the inverse map: \deqn{(f^{-1})'(z) = \exp(z)}}
#'   \item{invert_hessian(x, lower)}{Second derivative of the inverse map: \deqn{(f^{-1})''(z) = \exp(z)}}
#' }
#'
#' @return A list of six functions: \code{map}, \code{invert}, \code{map_jacobian},
#'   \code{map_hessian}, \code{invert_jacobian}, \code{invert_hessian}.
#'
#' @export
#'
#' @examples
#' f <- map_positive()
#' x <- 1:5
#' f$map(x, lower = 0)
#' f$invert(f$map(x, lower = 0), lower = 0)
#' f$map_jacobian(x, lower = 0)
#' f$invert_jacobian(f$map(x, lower = 0), lower = 0)
#' f$map_hessian(x, lower = 0)
#' f$invert_hessian(f$map(x, lower = 0), lower = 0)
map_positive <- function() {

  map <- function(x, lower) {
    log(x - lower)
  }

  invert <- function(x, lower) {
    lower + exp(x)
  }

  map_jacobian <- function(x, lower) {
    1/(x - lower)
  }

  map_hessian <- function(x, lower) {
    -1 / (x - lower)^2
  }

  invert_jacobian <- function(x, lower) {
    exp(x)
  }

  invert_hessian <- function(x, lower) {
    exp(x)
  }

  list(map = map, invert = invert,
       map_jacobian = map_jacobian, map_hessian = map_hessian,
       invert_jacobian = invert_jacobian, invert_hessian = invert_hessian)

}



#-------------------------------------------------------------------------------



#' Create mapping functions between half-open intervals and the real line
#'
#' Returns a list of functions implementing a bijective transformation between
#' a half-open interval `(-Inf, upper)` and the real line using a logarithmic transformation,
#' together with their derivatives and inverses.
#'
#' @details
#' The returned list contains:
#' \describe{
#'   \item{map(x, upper)}{Maps from (-Inf, upper) to the real line: \deqn{f(x) = \log(\texttt{upper} - x)}}
#'   \item{invert(x, upper)}{Inverse mapping: \deqn{f^{-1}(z) = \texttt{upper} - \exp(z)}}
#'   \item{map_jacobian(x, upper)}{First derivative of the forward map: \deqn{f'(x) = -1/(\texttt{upper} - x)}}
#'   \item{map_hessian(x, upper)}{Second derivative of the forward map: \deqn{f''(x) = -1/(\texttt{upper} - x)^2}}
#'   \item{invert_jacobian(x, upper)}{First derivative of the inverse map: \deqn{(f^{-1})'(z) = -\exp(z)}}
#'   \item{invert_hessian(x, upper)}{Second derivative of the inverse map: \deqn{(f^{-1})''(z) = -\exp(z)}}
#' }
#'
#' @return A list of six functions: \code{map}, \code{invert}, \code{map_jacobian},
#'   \code{map_hessian}, \code{invert_jacobian}, \code{invert_hessian}.
#'
#' @export
#'
#' @examples
#' f <- map_negative()
#' x <- -(5:1)
#' f$map(x, upper = 0)
#' f$invert(f$map(x, upper = 0), upper = 0)
#' f$map_jacobian(x, upper = 0)
#' f$invert_jacobian(f$map(x, upper = 0), upper = 0)
#' f$map_hessian(x, upper = 0)
#' f$invert_hessian(f$map(x, upper = 0), upper = 0)
map_negative <- function() {

  map <- function(x, upper) {
    log(upper - x)
  }

  invert <- function(x, upper) {
    upper - exp(x)
  }

  map_jacobian <- function(x, upper) {
    -1 / (upper - x)
  }

  map_hessian <- function(x, upper) {
    -1 / (upper - x)^2
  }

  invert_jacobian <- function(x, upper) {
    -exp(x)
  }

  invert_hessian <- function(x, upper) {
    -exp(x)
  }

  list(map = map, invert = invert,
       map_jacobian = map_jacobian, map_hessian = map_hessian,
       invert_jacobian = invert_jacobian, invert_hessian = invert_hessian)

}



#-------------------------------------------------------------------------------



#' Create mapping functions between constrained parameter spaces and the real line
#'
#' Returns a set of transformation functions that map parameters from various constrained spaces
#' to the real line and vice versa, along with their first and second derivatives.
#'
#' @param lower A numeric vector containing the lower bounds for each parameter.
#'   Use \code{-Inf} for unbounded below parameters.
#' @param upper A numeric vector containing the upper bounds for each parameter.
#'   Use \code{Inf} for unbounded above parameters.
#'
#' @details
#' The function automatically selects the appropriate transformation for each parameter based on its bounds:
#' \itemize{
#'   \item For parameters bounded both below and above \code{(lower[i], upper[i])}, [map_interval()] is used.
#'   \item For parameters bounded only below \code{(lower[i], Inf)}, [map_positive()] is used.
#'   \item For parameters bounded only above \code{(-Inf, upper[i])}, [map_negative()] is used.
#'   \item For completely unbounded parameters \code{(-Inf, Inf)}, the identity function is used.
#' }
#'
#' The returned list contains six functions:
#' \describe{
#'   \item{map(par)}{Maps parameters from the constrained space to the real line.}
#'   \item{invert(par)}{Maps parameters from the real line back to the constrained space.}
#'   \item{map_jacobian(par)}{Returns the first derivative (Jacobian) of the forward transformation.}
#'   \item{map_hessian(par)}{Returns the second derivative (Hessian) of the forward transformation.}
#'   \item{invert_jacobian(par)}{Returns the first derivative (Jacobian) of the inverse transformation.}
#'   \item{invert_hessian(par)}{Returns the second derivative (Hessian) of the inverse transformation.}
#' }
#' Each function expects a numeric vector \code{par} of the same length as \code{lower} and \code{upper}.
#'
#' @return A list containing six functions: \code{map}, \code{invert}, \code{map_jacobian},
#'   \code{map_hessian}, \code{invert_jacobian}, \code{invert_hessian}.
#'
#' @export
#'
#' @examples
#' # Define constraints for parameters
#' lower <- c(-Inf, 0, -1)
#' upper <- c(Inf, Inf, 1)
#'
#' # Create mapping functions
#' map_functions <- make_map_function(lower, upper)
#'
#' # Extract the individual functions
#' map <- map_functions$map
#' invert <- map_functions$invert
#' map_jacobian <- map_functions$map_jacobian
#' map_hessian <- map_functions$map_hessian
#' invert_jacobian <- map_functions$invert_jacobian
#' invert_hessian <- map_functions$invert_hessian
#'
#' # Define parameter values in the constrained space
#' x <- c(0, 3, 0.2)
#'
#' # Map to unconstrained space
#' y <- map(x)
#'
#' # Invert mapping back to constrained space
#' x_recovered <- invert(y)
#'
#' # Verify that the mapping is bijective
#' all.equal(x, x_recovered)
#'
#' # Compute derivatives
#' j_map <- map_jacobian(x)
#' h_map <- map_hessian(x)
#' j_invert <- invert_jacobian(y)
#' h_invert <- invert_hessian(y)
make_map_function <- function(lower, upper) {
  if (length(lower) != length(upper)) {
    stop("lower and upper must be of the same length")
  }

  if (any(lower >= upper)) stop("lower must be less than upper")

  npar <- length(lower)
  map_string <- character(npar)
  invert_string <- map_string
  map_j_string <- map_string
  map_h_string <- map_string
  invert_j_string <- map_string
  invert_h_string <- map_string

  for (i in 1:npar) {
    if (is.infinite(lower[i]) & is.infinite(upper[i])) {
      map_string[i] <- invert_string[i] <- paste0("par[", deparse(i), "]")
      map_j_string[i] <- invert_j_string[i] <- "1"
      map_h_string[i] <- invert_h_string[i] <- "0"
    } else if (!is.infinite(lower[i]) & is.infinite(upper[i])) {
      map_string[i] <- paste0("map_positive()$map(par[", deparse(i), "], ", deparse(lower[i]), ")")
      invert_string[i] <- paste0("map_positive()$invert(par[", deparse(i), "], ", deparse(lower[i]), ")")
      map_j_string[i] <- paste0("map_positive()$map_jacobian(par[", deparse(i), "], ", deparse(lower[i]), ")")
      map_h_string[i] <- paste0("map_positive()$map_hessian(par[", deparse(i), "], ", deparse(lower[i]), ")")
      invert_j_string[i] <- paste0("map_positive()$invert_jacobian(par[", deparse(i), "], ", deparse(lower[i]), ")")
      invert_h_string[i] <- paste0("map_positive()$invert_hessian(par[", deparse(i), "], ", deparse(lower[i]), ")")
    } else if (is.infinite(lower[i]) & !is.infinite(upper[i])) {
      map_string[i] <- paste0("map_negative()$map(par[", deparse(i), "], ", deparse(upper[i]), ")")
      invert_string[i] <- paste0("map_negative()$invert(par[", deparse(i), "], ", deparse(upper[i]), ")")
      map_j_string[i] <- paste0("map_negative()$map_jacobian(par[", deparse(i), "], ", deparse(upper[i]), ")")
      map_h_string[i] <- paste0("map_negative()$map_hessian(par[", deparse(i), "], ", deparse(upper[i]), ")")
      invert_j_string[i] <- paste0("map_negative()$invert_jacobian(par[", deparse(i), "], ", deparse(upper[i]), ")")
      invert_h_string[i] <- paste0("map_negative()$invert_hessian(par[", deparse(i), "], ", deparse(upper[i]), ")")
    } else if (!is.infinite(lower[i]) & !is.infinite(upper[i])) {
      map_string[i] <- paste0("map_interval()$map(par[", deparse(i), "], ", deparse(lower[i]), ", ", deparse(upper[i]), ")")
      invert_string[i] <- paste0("map_interval()$invert(par[", deparse(i), "], ", deparse(lower[i]), ", ", deparse(upper[i]), ")")
      map_j_string[i] <- paste0("map_interval()$map_jacobian(par[", deparse(i), "], ", deparse(lower[i]), ", ", deparse(upper[i]), ")")
      map_h_string[i] <- paste0("map_interval()$map_hessian(par[", deparse(i), "], ", deparse(lower[i]), ", ", deparse(upper[i]), ")")
      invert_j_string[i] <- paste0("map_interval()$invert_jacobian(par[", deparse(i), "], ", deparse(lower[i]), ", ", deparse(upper[i]), ")")
      invert_h_string[i] <- paste0("map_interval()$invert_hessian(par[", deparse(i), "], ", deparse(lower[i]), ", ", deparse(upper[i]), ")")
    }
  }

  map <- function(par) {
    code <- paste0("c(", paste(map_string, collapse = ", "), ")")
    eval(parse(text = code))
  }

  invert <- function(par) {
    code <- paste0("c(", paste(invert_string, collapse = ", "), ")")
    eval(parse(text = code))
  }

  map_jacobian <- function(par) {
    code <- paste0("c(", paste(map_j_string, collapse = ", "), ")")
    eval(parse(text = code))
  }

  map_hessian <- function(par) {
    code <- paste0("c(", paste(map_h_string, collapse = ", "), ")")
    eval(parse(text = code))
  }

  invert_jacobian <- function(par) {
    code <- paste0("c(", paste(invert_j_string, collapse = ", "), ")")
    eval(parse(text = code))
  }

  invert_hessian <- function(par) {
    code <- paste0("c(", paste(invert_h_string, collapse = ", "), ")")
    eval(parse(text = code))
  }

  list(
    map = map, invert = invert, map_jacobian = map_jacobian,
    map_hessian = map_hessian, invert_jacobian = invert_jacobian,
    invert_hessian = invert_hessian
  )
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
#' map_jacobian <- map_functions$map_jacobian
#'
#' # Map parameters to unbounded space
#' mapped_pars <- map(par)
#'
#' # Calculate Jacobian adjustment for covariance matrix based on delta method
#' j <- diag(map_jacobian(par))
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
  # if(!isSymmetric(v)) stop("v must be a symmetric matrix")
  if(nrow(v) != npar) stop("dimension of variance-covariance matrix and par must coincide")
  if(length(mult) != 1 & length(mult) != npar) stop("mult must be of dimension 1 or npar")
  eps <- tryCatch(c(mvtnorm::rmvnorm(1, sigma = v)),
                  error = function(e) {
                    stats::rnorm(npar, mean = 0, sd = sqrt(diag(v)))
                  }
  )
  par + mult*eps
}














