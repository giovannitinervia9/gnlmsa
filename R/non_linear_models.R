


#-------------------------------------------------------------------------------



#' Create a Linear Function and its Derivatives
#'
#' This function creates a linear function and its derivatives (Jacobian and Hessian) for use in optimization, modeling, and statistical analysis.
#'
#' @details
#' The linear function implemented here has the following form:
#' \deqn{f(x_i, \theta) = x_i^{\intercal}\theta}
#'
#' Where:
#' \itemize{
#'   \item \eqn{x_i^{\intercal}} is the \eqn{i}-th row vector of the matrix of covariates \eqn{X}, with each row representing an observation
#'         and each column representing a covariate
#'   \item \eqn{\theta} is a vector of coefficients
#'   \item The function returns a vector of linear combinations, one for each row of \eqn{X}
#' }
#'
#' The Jacobian (first derivatives) matrix has the form:
#' \deqn{J_{ij} = \dfrac{\partial y_i}{\partial \theta_j} =  x_{ij}}
#'
#' Where:
#' \itemize{
#'   \item \eqn{x_{ij}} is the value of covariate \eqn{j} for observation \eqn{i}
#' }
#'
#' The Hessian (second derivatives) is a zero matrix of dimension \eqn{k \times k} where \eqn{k} is the number of columns of \eqn{X}.
#'
#' Additionally, the function allows the inclusion of an intercept term and constraints on the coefficients (lower and upper bounds).
#'
#' @param X A matrix or data frame of covariates (input data). Each row represents an observation and each column represents a covariate.
#' @param intercept A logical indicating whether to include an intercept term (default is TRUE).
#' @param lower A vector of lower bounds for the coefficients. If not specified, defaults to \code{-Inf} for all coefficients.
#' @param upper A vector of upper bounds for the coefficients. If not specified, defaults to \code{Inf} for all coefficients.
#'
#' @return A list containing three functions and additional information:
#'   \item{f}{The linear function. Takes parameters:
#'     \itemize{
#'       \item \code{X}: A matrix or data frame of covariates
#'       \item \code{theta}: A vector of coefficients
#'     }
#'   }
#'   \item{J}{The Jacobian (first derivatives) function. Takes parameters:
#'     \itemize{
#'       \item \code{X}: A matrix or data frame of covariates
#'       \item \code{theta}: A vector of coefficients
#'       \item \code{y}: Optional; pre-computed output values from the linear function
#'     }
#'   }
#'   \item{H}{The Hessian (second derivatives) function. Takes parameters:
#'     \itemize{
#'       \item \code{X}: A matrix or data frame of covariates
#'       \item \code{theta}: A vector of coefficients
#'       \item \code{y}: Optional; pre-computed output values from the linear function
#'     }
#'     Returns a list of Hessian matrices, one for each observation in \code{X}
#'   }
#'   \item{lower}{A vector of lower bounds for the coefficients (if specified).}
#'   \item{upper}{A vector of upper bounds for the coefficients (if specified).}
#'   \item{X}{The input matrix or data frame of covariates.}
#'
#' @export
Linear <- function(X, intercept = T, lower, upper){
  f <- function(X, theta) drop(as.matrix(X)%*%theta)

  J <- function(X, theta){
    as.matrix(X)
  }

  H <- function(X, theta){
    k <- length(theta)
    hi <- matrix(0, k, k)
    rep(list(hi), NROW(X))
  }

  if (missing(X)) {
    return(list(f = f, J = J, H = H))
  }

  if (intercept) {
    if(!all(X[, 1] == 1)) X <- cbind(1, X)
  }

  if (missing(lower)) lower <- rep(-Inf, ncol(X))
  if (missing(upper)) upper <- rep(Inf, ncol(X))

  list(f = f, J = J, H = H, lower = lower, upper = upper, X = X)

}



#-------------------------------------------------------------------------------



#' Create an Exponential Function and its Derivatives
#'
#' This function creates an exponential function and its derivatives (Jacobian and Hessian) for use in optimization, modeling, and statistical analysis.
#'
#' @details
#' The exponential function implemented here has the following form:
#' \deqn{f(x_i, \theta) = \exp(x_i^\intercal \theta)}
#'
#' Where:
#' \itemize{
#'   \item \eqn{x_i^{\intercal}} is the \eqn{i}-th row vector of the matrix of covariates \eqn{X}, with each row representing an observation
#'         and each column representing a covariate
#'   \item \eqn{\theta} is a vector of coefficients
#'   \item The function returns a vector of exponential values, one for each row of \eqn{X}
#' }
#'
#' The Jacobian matrix has the form:
#' \deqn{J_{ij} = \dfrac{\partial y_i}{\partial \theta_j} =  y_i x_{ij}}
#'
#' Where:
#' \itemize{
#'   \item \eqn{y_i} is the output from the exponential function for observation \eqn{i}
#'   \item \eqn{x_{ij}} is the value of covariate \eqn{j} for observation \eqn{i}
#' }
#'
#' The Hessian matrix for the generic \eqn{i}-th observation has the structure:
#' \deqn{H_{jl} = H_{lj} = \dfrac{\partial^2y_i}{\partial \theta_j \partial \theta_l} = y_i x_{ij} x_{il}}
#'
#' Where:
#' \itemize{
#'   \item \eqn{y_i} is the output from the exponential function for the \eqn{i}-th observation
#'   \item \eqn{x_{ij}} and \eqn{x_{il}} are values of covariates \eqn{j} and \eqn{l} for the \eqn{i}-th observation
#' }
#'
#' Additionally, the function allows the inclusion of an intercept term and constraints on the coefficients (lower and upper bounds).
#'
#' @param X A matrix or data frame of covariates (input data). Each row represents an observation and each column represents a covariate.
#' @param intercept A logical indicating whether to include an intercept term (default is TRUE).
#' @param lower A vector of lower bounds for the coefficients. If not specified, defaults to \code{-Inf} for all coefficients.
#' @param upper A vector of upper bounds for the coefficients. If not specified, defaults to \code{Inf} for all coefficients.
#'
#' @return A list containing three functions and additional information:
#'   \item{f}{The exponential function. Takes parameters:
#'     \itemize{
#'       \item \code{X}: A matrix or data frame of covariates
#'       \item \code{theta}: A vector of coefficients
#'     }
#'   }
#'   \item{J}{The Jacobian (first derivatives) function. Takes parameters:
#'     \itemize{
#'       \item \code{X}: A matrix or data frame of covariates
#'       \item \code{theta}: A vector of coefficients
#'       \item \code{y}: Optional; pre-computed output values from the exponential function
#'     }
#'   }
#'   \item{H}{The Hessian (second derivatives) function. Takes parameters:
#'     \itemize{
#'       \item \code{X}: A matrix or data frame of covariates
#'       \item \code{theta}: A vector of coefficients
#'       \item \code{y}: Optional; pre-computed output values from the exponential function
#'     }
#'     Returns a list of Hessian matrices, one for each observation in \code{X}
#'   }
#'   \item{lower}{A vector of lower bounds for the coefficients (if specified).}
#'   \item{upper}{A vector of upper bounds for the coefficients (if specified).}
#'   \item{X}{The input matrix or data frame of covariates.}
#'
#' @export
Exp <- function(X, intercept = TRUE, lower, upper) {
  f <- function(X, theta) {
    drop(exp(as.matrix(X) %*% theta))
  }

  J <- function(X, theta, y = NULL) {
    X <- as.matrix(X)
    if (is.null(y)) {
      y <- exp(X %*% theta)
    }
    X * as.vector(y)
  }

  H <- function(X, theta, y = NULL) {
    X <- as.matrix(X)
    if (is.null(y)) {
      y <- exp(X %*% theta)
    }

    k <- length(theta)
    n <- nrow(X)

    # Preallocate the list of Hessians
    H_list <- vector("list", n)

    # Efficient row-wise Hessian computation (no Map or split)
    for (i in seq_len(n)) {
      xi <- X[i, ]
      yi <- y[i]
      H_list[[i]] <- tcrossprod(xi) * yi
    }

    H_list
  }

  if (missing(X)) {
    return(list(f = f, J = J, H = H))
  }

  if (intercept) {
    if(!all(X[, 1] == 1)) X <- cbind(1, X)
  }

  if (missing(lower)) lower <- rep(-Inf, ncol(X))
  if (missing(upper)) upper <- rep(Inf, ncol(X))

  list(f = f, J = J, H = H, lower = lower, upper = upper, X = X)
}



#-------------------------------------------------------------------------------



#' Create a Cobb-Douglas Production Function and its Derivatives
#'
#' This function creates a Cobb-Douglas production function and its derivatives
#' (Jacobian and Hessian) for use in optimization and modeling. The Cobb-Douglas
#' production function is commonly used in economics to represent the relationship
#' between inputs and output levels.
#'
#' @details
#' The Cobb-Douglas production function has the following general form:
#' \deqn{f(x_i, \theta) = \theta_0 \prod_{j=1}^k x_{ij}^{\intercal \theta_{j}}}
#'
#' Where:
#' \itemize{
#'   \item \eqn{x_i^{\intercal}} is the \eqn{i}-th row vector of the matrix of input factors \eqn{X}, with each row representing an observation
#'         and each column representing an input factor
#'   \item \eqn{k} is the number of input factors
#'   \item \eqn{\theta_0} is a scale parameter (efficiency parameter)
#'   \item \eqn{\theta_{j}} are the elasticity parameters for each input factor
#'   \item The function exhibits constant returns to scale when \eqn{\sum_{j=1}^k \theta_j = 1}
#' }
#'
#' The Jacobian (first derivatives) matrix has the form:
#' \deqn{J_{i1} = \dfrac{\partial y_i}{\partial \theta_0} = \dfrac{y_i}{\theta_0}}
#' \deqn{J_{ij} = \dfrac{\partial y_i}{\partial \theta_j} =  y_i \log(x_{ij}) \text{ for } j = 1, \ldots, k}
#'
#' Where:
#' \itemize{
#'   \item \eqn{x_{ij}} is the value of input factor \eqn{j} for the \eqn{i}-th observation
#'   \item \eqn{y} is the output from the Cobb-Douglas function
#'   \item \eqn{k} is the number of input factors
#' }
#'
#' The Hessian (second derivatives) matrix for the \eqn{i}-th observation has the form:
#' \deqn{H_{11} = \dfrac{\partial^2 y_i}{\partial ^2 \theta_0} = 0}
#' \deqn{H_{1j} = H_{j,1} = \dfrac{\partial ^2y_i}{\partial \theta_0 \partial \theta_j} = \dfrac{y_i}{\theta_0} \log(x_{ij}) \text{ for } j = 1, \ldots, k}
#' \deqn{H_{jl} = H_{l,j} = \dfrac{\partial ^2 y_i}{\partial \theta_j \partial \theta_l} = y_i \log(x_{ij}) \log(x_{il}) \text{ for } j,l = 1, \ldots, k}
#'
#' Additionally, the function allows the inclusion of constraints on the coefficients (lower and upper bounds).
#'
#' @param X A matrix or data frame of input factors (input data). Each row represents an observation and each column represents an input factor.
#' @param lower A vector of lower bounds for the coefficients. If not specified, defaults to \code{0} for the parameters.
#' @param upper A vector of upper bounds for the coefficients. If not specified, defaults to \code{Inf} for the parameters.
#'
#' @return A list containing three functions and additional information:
#'   \item{f}{The Cobb-Douglas production function. Takes parameters:
#'     \itemize{
#'       \item \code{X}: A matrix or data frame of input factors
#'       \item \code{theta}: A vector of parameters where \code{theta[1]} is the efficiency parameter and
#'             \code{theta[2:length(theta)]} are the elasticity parameters for each input factor
#'     }
#'   }
#'   \item{J}{The Jacobian (first derivatives) function. Takes parameters:
#'     \itemize{
#'       \item \code{X}: A matrix or data frame of input factors
#'       \item \code{theta}: A vector of parameters
#'       \item \code{y}: Optional; pre-computed output values from the Cobb-Douglas function
#'     }
#'   }
#'   \item{H}{The Hessian (second derivatives) function. Takes parameters:
#'     \itemize{
#'       \item \code{X}: A matrix or data frame of input factors
#'       \item \code{theta}: A vector of parameters
#'       \item \code{y}: Optional; pre-computed output values from the Cobb-Douglas function
#'     }
#'     Returns a list of Hessian matrices, one for each observation in \code{X}
#'   }
#'   \item{lower}{A vector of lower bounds for the coefficients (if specified).}
#'   \item{upper}{A vector of upper bounds for the coefficients (if specified).}
#'   \item{X}{The input matrix or data frame of input factors.}#'
#' @export
cobb_douglas <- function(X, lower, upper) {

  f <- function(X, theta) {
    drop(theta[1] * exp(log(as.matrix(X)) %*% theta[-1]))
  }

  J <- function(X, theta, y = NULL) {
    logX <- log(as.matrix(X))
    if (is.null(y)) {
      y <- drop(theta[1] * exp(logX %*% theta[-1]))
    }

    cbind(y / theta[1], y * logX)
  }

  H <- function(X, theta, y = NULL) {
    logX <- log(as.matrix(X))
    k <- length(theta)

    if (is.null(y)) {
      y <- drop(theta[1] * exp(logX %*% theta[-1]))
    }

    n <- nrow(X)
    H_list <- vector("list", n)

    for (i in seq_len(n)) {
      h <- matrix(0, k, k)
      lxi <- drop(logX[i, ])
      yi <- y[i]

      # Derivatives wrt theta_0 (intercept) and rest
      h[1, -1] <- h[-1, 1] <- (yi / theta[1]) * lxi
      for (r in 2:k) {
        for (c in r:k) {
          h[r, c] <- h[c, r] <- yi * lxi[r - 1] * lxi[c - 1]
        }
      }

      H_list[[i]] <- h
    }

    H_list
  }

  if (missing(X)) {
    return(list(f = f, J = J, H = H))
  }

  if (missing(lower)) lower <- rep(0, ncol(X) + 1)
  if (missing(upper)) upper <- rep(Inf, ncol(X) + 1)
  list(f = f, J = J, H = H, lower = lower, upper = upper, X = X)
}


#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
