


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
#' @return A list containing three functions:
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
#'
#'
#' @examples
#' # Create the exponential function
#' lin_func <- Linear()
#'
#' # Define covariate data
#' X <- matrix(c(1, 0.5, 0.2,
#'               1, 0.8, 0.3,
#'               1, 0.6, 0.4), ncol = 3, byrow = TRUE)
#'
#' # Define coefficients
#' theta <- c(0.1, 0.5, -0.3)
#'
#' # Calculate linear combinations
#' y_values <- lin_func$f(X, theta)
#' y_values
#'
#' # Calculate Jacobian (first derivatives)
#' jacobian <- lin_func$J(X, theta)
#' jacobian
#'
#' # Calculate Hessian matrices (second derivatives)
#' hessian_list <- lin_func$H(X, theta)
#' hessian_list[[1]]  # Hessian matrix for the first observation
#'
#' @export
Linear <- function(){
  f <- function(X, theta) drop(X%*%theta)

  J <- function(X, theta, y){
    if(missing(y)) y <- Linear()$f(X, theta)
    X
  }

  H <- function(X, theta){
    k <- length(theta)
    hi <- matrix(0, k, k)
    rep(list(hi), NROW(X))
  }

  list(f = f, J = J, H = H)

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
#' The Hessian matrix for the generic \eqn{i}-th observation has the sctructure:
#' \deqn{H_{jl} = H_{lj} = \dfrac{\partial^2y_i}{\partial \theta_j \partial \theta_l} = y_i x_{ij} x_{il}}
#'
#' Where:
#' \itemize{
#'   \item \eqn{y_i} is the output from the exponential function for the \eqn{i}-th observation
#'   \item \eqn{x_{ij}} and \eqn{x_{il}} are values of covariates \eqn{j} and \eqn{l} for the \eqn{i}-th observation
#' }
#'
#' @return A list containing three functions:
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
#'
#'
#' @examples
#' # Create the exponential function
#' exp_func <- Exp()
#'
#' # Define covariate data
#' X <- matrix(c(1, 0.5, 0.2,
#'               1, 0.8, 0.3,
#'               1, 0.6, 0.4), ncol = 3, byrow = TRUE)
#'
#' # Define coefficients
#' theta <- c(0.1, 0.5, -0.3)
#'
#' # Calculate exponential values
#' y_values <- exp_func$f(X, theta)
#' y_values
#'
#' # Calculate Jacobian (first derivatives)
#' jacobian <- exp_func$J(X, theta)
#' jacobian
#'
#' # Calculate Hessian matrices (second derivatives)
#' hessian_list <- exp_func$H(X, theta)
#' hessian_list[[1]]  # Hessian matrix for the first observation
#'
#' @export
Exp <- function(){
  f <- function(X, theta) drop(exp(X%*%theta))

  J <- function(X, theta, y) {
    if (missing(y)) {
      y <- Exp()$f(X, theta)
    }
    y*X
  }

  H <- function(X, theta, y) {
    if (missing(y)) {
      y <- Exp()$f(X, theta)
    }
    k <- length(theta)

    hi <- function(x, y) {
      h <- matrix(NA, k, k)
      for(i in 1:k){
        for(j in i:k){
          h[i, j] <- h[j, i] <- y*x[i]*x[j]
        }
      }
      h
    }

    Map(function(x, y) hi(x, y), x = split(as.matrix(X), f = row(as.matrix(X))), y = y)
  }

  list(f = f, J = J, H = H)

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
#'         and each column representing a input factor
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
#' @return A list containing three functions:
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
#'
#' @examples
#' # Create the Cobb-Douglas function
#' cd <- cobb_douglas()
#'
#' # Define input data (labor and capital)
#' X <- data.frame(labor = c(10, 20, 15), capital = c(50, 40, 60))
#'
#' # Define parameters: efficiency and elasticities
#' # theta[1] = efficiency, theta[2] = labor elasticity, theta[3] = capital elasticity
#' theta <- c(2, 0.6, 0.4)  # Constant returns to scale as 0.6 + 0.4 = 1
#'
#' # Calculate production output
#' output <- cd$f(X, theta)
#' output
#'
#' # Calculate Jacobian (first derivatives)
#' jacobian <- cd$J(X, theta)
#' jacobian
#'
#' # Calculate Hessian matrices (second derivatives)
#' hessian_list <- cd$H(X, theta)
#' hessian_list[[1]]  # Hessian matrix for the first observation
#'
#' @export
cobb_douglas <- function(){

  f <- function(X, theta) {
    apply(as.matrix(X), 1, function(x) theta[1]*prod(x^theta[-1]))
  }


  J <- function(X, theta, y) {
    if (missing(y)) {
      y <- apply(as.matrix(X), 1, function(x) theta[1]*prod(x^theta[-1]))
    }

    j <- cbind(y/theta[1])

    for(i in 1:ncol(X)){
      j <- cbind(j, y*log(X[, i]))
    }
    j
  }


  H <- function(X, theta, y){
    if (missing(y)) {
      y <- apply(as.matrix(X), 1, function(x) theta[1]*prod(x^theta[-1]))
    }
    k <- length(theta)

    hi <- function(lx, y) {
      h <- matrix(NA, k, k)
      for(i in 1:k){
        for(j in i:k){

          if ((i == 1) & (j == 1)) {
            h[i, j]  <- 0
          }
          else if (i == 1) {
            h[i, j] <- h[j, i] <- (y/theta[1])*lx[j - 1]
          }
          else {
            h[i, j] <- h[j, i] <-  y*lx[i - 1]*lx[j - 1]
          }
        }

      }

      h
    }

    Map(function(x, y) hi(x, y), x = split(log(as.matrix(X)), f = row(as.matrix(X))), y = y)

  }

  list(f = f, J = J, H = H)


}



#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
