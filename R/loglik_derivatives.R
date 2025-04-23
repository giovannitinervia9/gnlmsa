


#-------------------------------------------------------------------------------



#' Gradient with Respect to Mean Component Parameters in Generalized Non-Linear Models
#'
#' Computes the gradient of the log-likelihood with respect to the parameters of the mean component
#' in a Generalized Non-Linear Model (GNLM), assuming the response variable belongs to an exponential family.
#'
#' @param y Numeric vector of observed responses.
#' @param X Numeric matrix of covariates for the mean component.
#' @param beta Numeric vector of current values for the parameters of the mean component.
#' @param mu Numeric vector of conditional means, typically obtained via the inverse link function applied to \code{eta}.
#' @param eta Numeric vector of nonlinear predictors \eqn{\eta = \eta(X, \beta)}.
#' @param phi Numeric vector of conditional dispersion values.
#' @param f_mu A user-defined function of signature \code{f(X, theta)} returning the nonlinear predictor.
#' @param J_mu Optional. A function to compute the Jacobian of \code{f_mu}, with signature \code{J_mu(X, theta)}.
#'   Returns a \eqn{n \times k} matrix, where `n = length(y)` and `k = length(beta)`.
#' @param mu.eta Function computing the derivative \eqn{\partial \mu / \partial \eta}, as returned by [make_link()].
#' @param variance Function defining the variance of the distribution as a function of \code{mu} and \code{phi}.
#'
#' @return A numeric vector of length equal to the number of parameters in the mean component.
#'
#' @details
#' This function evaluates the score vector (first derivative of the log-likelihood) with respect to
#' the mean parameters \eqn{\beta}. It is valid for any GNLM family provided that the appropriate
#' variance and derivative functions are specified.
#'
#'
#'
#' @examples
#' \dontrun{
#' # create family object
#' fam_gamma <- gnlmsa_Gamma(link_mu = "identity", link_phi = "log")
#'
#' # create non linear function for mean
#' cd <- cobb_douglas()
#' f_mu <- cd$f
#' J_mu <- cd$J
#'
#' # create non linear function for dispersion
#' lin <- Linear()
#' f_phi <- lin$f
#' J_phi <- lin$J
#'
#' # sample size
#' n <- 100
#'
#' # covariates for the mean component
#' X <- cbind(rgamma(n, 1000, 10), rgamma(n, 2000, 5))
#'
#' # covariates for the variance component
#' Z <- cbind(1, X)
#'
#' # parameters for mean component
#' beta <- c(1, .5, .3)
#' eta <- f_mu(X, beta)
#' mu <- fam_gamma$linkinv_mu(eta)
#'
#' # parameters for dispersion component
#' gamma <- c(1, .2, -.05)
#' vi <- f_phi(Z, gamma)
#' phi <- fam_gamma$linkinv_phi(vi)
#'
#' # simulate response
#' y <- rgamma(n, phi, phi/mu)
#'
#' grad_mu(y, X, beta, mu, eta, phi,
#'         f_mu, J_mu,
#'         fam_gamma$mu.eta, fam_gamma$variance)
#' }
#' @export
grad_mu <- function(y, X, beta, mu, eta, phi, f_mu, J_mu, mu.eta, variance){
  if(missing(J_mu)) J_mu <- make_jacobian(f_mu)
  colSums(J_mu(X, beta)*(mu.eta(eta)/variance(mu, phi))*(y - mu)) # J*w*(y - mu)
}



#-------------------------------------------------------------------------------



#' Hessian with Respect to Mean Component Parameters in Generalized Non-Linear Models
#'
#' Computes the Hessian matrix (second derivative of the log-likelihood) with respect to the parameters
#' of the mean component in a Generalized Non-Linear Model (GNLM), assuming an exponential family distribution.
#'
#' @param y Numeric vector of observed responses.
#' @param X Numeric matrix of covariates associated with the mean component.
#' @param beta Numeric vector of current values for the mean component parameters.
#' @param mu Numeric vector of conditional means, typically computed as the inverse link applied to \code{eta}.
#' @param eta Numeric vector of nonlinear predictors \eqn{\eta = \eta(X, \beta)}.
#' @param phi Numeric vector of conditional dispersion values.
#' @param f_mu Function with signature \code{f(X, theta)} that computes the nonlinear predictor for the mean component.
#' @param J_mu Optional. Function computing the Jacobian of \code{f_mu}; must return an \eqn{n \times k} matrix.
#' @param H_mu Optional. Function computing the Hessian of \code{f_mu}; must return a list of \eqn{n} Hessian matrices,
#'   one per observation.
#' @param mu.eta Function computing \eqn{\partial \mu / \partial \eta}, as returned by [make_link()].
#' @param mu2.eta2 Function computing \eqn{\partial^2 \mu / \partial^2 \eta^2}, as returned by [make_link()].
#' @param variance Function returning the conditional variance of the response as a function of \code{mu} and \code{phi}.
#' @param expected Logical. If \code{TRUE} (default), returns the expected Hessian.
#'   If \code{FALSE}, returns the observed (empirical) Hessian.
#'
#' @return A square symmetric matrix of dimension equal to the number of parameters in the mean component.
#'
#' @details
#' This function evaluates the hessian matrix (second derivative of the log-likelihood) with respect to
#' the mean parameters \eqn{\beta}. It is valid for any GNLM family provided that the appropriate
#' variance and derivative functions are specified.
#' @export
#' @examples
#' # create family object
#' fam_gamma <- gnlmsa_Gamma(link_mu = "identity", link_phi = "log")
#'
#' # create non linear function for mean
#' cd <- cobb_douglas()
#' f_mu <- cd$f
#' J_mu <- cd$J
#' H_mu <- cd$H
#'
#' # create non linear function for dispersion
#' lin <- Linear()
#' f_phi <- lin$f
#'
#' # sample size
#' n <- 100
#'
#' # covariates for the mean component
#' X <- cbind(rgamma(n, 1000, 10), rgamma(n, 2000, 5))
#'
#' # covariates for the variance component
#' Z <- cbind(1, X)
#'
#' # parameters for mean component
#' beta <- c(1, .5, .3)
#' eta <- f_mu(X, beta)
#' mu <- fam_gamma$linkinv_mu(eta)
#'
#' # parameters for dispersion component
#' gamma <- c(1, .2, -.05)
#' vi <- f_phi(Z, gamma)
#' phi <- fam_gamma$linkinv_phi(vi)
#'
#' # simulate response
#' y <- rgamma(n, phi, phi/mu)
#'
#' hess_mu(y, X, beta, mu, eta, phi,
#'         f_mu, J_mu, H_mu,
#'         fam_gamma$mu.eta, fam_gamma$mu2.eta2, fam_gamma$variance,
#'         expected = TRUE)
#' hess_mu(y, X, beta, mu, eta, phi,
#'         f_mu, J_mu, H_mu,
#'         fam_gamma$mu.eta, fam_gamma$mu2.eta2, fam_gamma$variance,
#'         expected = FALSE)
hess_mu<- function(y, X, beta, mu, eta, phi,
                   f_mu, J_mu, H_mu,
                   mu.eta, mu2.eta2, variance,
                   expected = TRUE) {

  if (missing(J_mu)) J_mu <- make_jacobian(f_mu)
  if (missing(H_mu)) H_mu <- make_hessian(f_mu)

  mu_eta <- mu.eta(eta)
  mu2_eta2 <- mu2.eta2(eta)

  j <- J_mu(X, beta)
  v <- variance(mu, phi)

  # Expected part
  w1 <- (mu_eta^2) / v
  h1 <- -crossprod(j, w1 * j)

  if (expected) {
    return(h1)
  }

  # Observed part
  w2 <- (y - mu) / v

  h21 <- crossprod(j, (w2 * mu2_eta2) * j)

  H_list <- H_mu(X, beta)
  h22 <- matrix(0, ncol(j), ncol(j))
  for (i in seq_along(w2)) {
    h22 <- h22 + w2[i] * mu_eta[i] * H_list[[i]]
  }

  h1 + h21 + h22
}



#-------------------------------------------------------------------------------



#' Construct the Full Hessian Matrix for a Generalized Non-Linear Model
#'
#' Assembles the full Hessian matrix of the log-likelihood function for a Generalized Non-Linear Model (GNLM),
#' by combining the individual Hessians of the mean and dispersion components along with their cross-derivatives.
#'
#' @param h_mu Hessian matrix of the log-likelihood with respect to the parameters of the mean component.
#' @param h_phi Hessian matrix of the log-likelihood with respect to the parameters of the dispersion component.
#' @param h_mu_phi Cross-Hessian matrix of the log-likelihood with respect to mean and dispersion parameters.
#'
#' @return A symmetric square matrix of dimension \eqn{(k + p) \times (k + p)}, where:
#' \itemize{
#'   \item \eqn{k} is the number of parameters in the mean component,
#'   \item \eqn{p} is the number of parameters in the dispersion component.
#' }
#' The resulting matrix has the block structure:
#' \deqn{
#' \begin{bmatrix}
#' H_{\beta \beta} & H_{\beta \gamma} \\
#' H_{\beta \gamma}^\intercal & H_{\gamma \gamma}
#' \end{bmatrix}
#' }
#'
#' @details
#' This function is useful for second-order optimization methods (e.g., Newton-Raphson),
#' model diagnostics, or computing standard errors when both components of the model
#' (mean and dispersion) are estimated jointly.
#'
#' All inputs must be compatible in dimension and properly computed; no internal checks are performed.
#' @export
#' @examples
#' # create family object
#' fam_gamma <- gnlmsa_Gamma(link_mu = "identity", link_phi = "log")
#'
#' # create non linear function for mean
#' cd <- cobb_douglas()
#' f_mu <- cd$f
#' J_mu <- cd$J
#' H_mu <- cd$H
#'
#' # create non linear function for dispersion
#' lin <- Linear()
#' f_phi <- lin$f
#' J_phi <- lin$J
#' H_phi <- lin$H
#'
#' # sample size
#' n <- 100
#'
#' # covariates for the mean component
#' X <- cbind(rgamma(n, 1000, 10), rgamma(n, 2000, 5))
#'
#' # covariates for the variance component
#' Z <- cbind(1, X)
#'
#' # parameters for mean component
#' beta <- c(1, .5, .3)
#' eta <- f_mu(X, beta)
#' mu <- fam_gamma$linkinv_mu(eta)
#'
#' # parameters for dispersion component
#' gamma <- c(1, .2, -.05)
#' vi <- f_phi(Z, gamma)
#' phi <- fam_gamma$linkinv_phi(vi)
#'
#' # simulate response
#' y <- rgamma(n, phi, phi/mu)
#'
#' # compute h_mu
#' h_mu <- hess_mu(y, X, beta, mu, eta, phi,
#'                 f_mu, J_mu, H_mu,
#'                 fam_gamma$mu.eta, fam_gamma$mu2.eta2, fam_gamma$variance)
#'
#' # compute h_phi
#' h_phi <- fam_gamma$hess_phi(y, Z, gamma,
#'                             phi, vi, mu,
#'                             f_phi, J_phi, H_phi,
#'                             fam_gamma$phi.vi, fam_gamma$phi2.vi2)
#'
#' # compute h_mu_phi
#' h_mu_phi <- fam_gamma$hess_mu_phi(y, X, Z, beta, gamma,
#'                                   mu, eta, phi, vi,
#'                                   f_mu, J_mu, f_phi, J_phi,
#'                                   fam_gamma$mu.eta, fam_gamma$phi.vi)
#'
#' # build full hessian
#' hess(h_mu, h_phi, h_mu_phi)
hess <- function(h_mu, h_phi, h_mu_phi){
  rbind(cbind(h_mu, h_mu_phi), cbind(t(h_mu_phi), h_phi))
}



#-------------------------------------------------------------------------------



#' Compute the Negative Log-Likelihood for a Generalized Non-Linear Model
#'
#' Evaluates the negative log-likelihood function for a Generalized Non-Linear Model (GNLM)
#' given a parameter vector and a user-supplied log-likelihood function from a `family_gnlmsa` object.
#' This function is typically used as an objective function in optimization routines.
#'
#' @param y Numeric vector of observed responses.
#' @param X Numeric matrix of covariates for the mean component.
#' @param Z Numeric matrix of covariates for the dispersion component.
#' @param par Numeric vector containing all model parameters. The first \code{npar_mu} elements are
#'   assumed to be associated with the mean component, and the remainder with the dispersion component.
#' @param npar_mu Integer. Number of parameters associated with the mean component.
#' @param loglik Function that computes the observation-wise log-likelihood contributions, as provided by a `gnlmsa_*()` family constructor.
#' @param f_mu Function with signature \code{f(X, theta)} computing the nonlinear predictor for the mean component.
#' @param f_phi Function with signature \code{f(Z, gamma)} computing the nonlinear predictor for the dispersion component.
#' @param linkinv_mu Inverse link function for the mean component (e.g., from [make_link()]).
#' @param linkinv_phi Inverse link function for the dispersion component.
#'
#' @return A single numeric value: the negative log-likelihood evaluated at the parameter vector \code{par}.
#'
#' @details
#' This function separates the full parameter vector \code{par} into \code{beta} (for the mean)
#' and \code{gamma} (for the dispersion), computes the predicted mean and dispersion using the
#' inverse link functions, and then evaluates the log-likelihood. The log-likelihood is returned with a
#' negative sign for compatibility with minimization algorithms.
#' @export
#' @examples
#' # create family object
#' fam_gamma <- gnlmsa_Gamma(link_mu = "identity", link_phi = "log")
#'
#' # create non linear function for mean
#' cd <- cobb_douglas()
#' f_mu <- cd$f
#' J_mu <- cd$J
#' H_mu <- cd$H
#'
#' # create non linear function for dispersion
#' lin <- Linear()
#' f_phi <- lin$f
#'
#' # sample size
#' n <- 100
#'
#' # covariates for the mean component
#' X <- cbind(rgamma(n, 1000, 10), rgamma(n, 2000, 5))
#'
#' # covariates for the variance component
#' Z <- cbind(1, X)
#'
#' # parameters for mean component
#' beta <- c(1, .5, .3)
#' eta <- f_mu(X, beta)
#' mu <- fam_gamma$linkinv_mu(eta)
#'
#' # parameters for dispersion component
#' gamma <- c(1, .2, -.05)
#' vi <- f_phi(Z, gamma)
#' phi <- fam_gamma$linkinv_phi(vi)
#'
#' # simulate response
#' y <- rgamma(n, phi, phi/mu)
#'
#' # compute log_loss
#' par <- c(beta, gamma)
#' npar_mu <- length(beta)
#' log_loss(y, X, Z, par, npar_mu, fam_gamma$loglik,
#'          f_mu, f_phi, fam_gamma$linkinv_mu, fam_gamma$linkinv_phi)
log_loss <- function(y, X, Z, par, npar_mu, loglik, f_mu, f_phi, linkinv_mu, linkinv_phi){
  beta <- par[1:npar_mu]
  gamma <- par[(npar_mu + 1):length(par)]
  mu <- linkinv_mu(f_mu(X, beta))
  phi <- linkinv_phi(f_phi(Z, gamma))
  -sum(loglik(y, mu, phi))
}



#-------------------------------------------------------------------------------



