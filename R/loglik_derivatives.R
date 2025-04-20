


#-------------------------------------------------------------------------------



#' Compute gradient with respect to parameters of mean component for Generalized Non-Linear Models
#'
#' @param y numerical vector of response variable.
#' @param X matrix containing the covariates for the mean component.
#' @param beta parameters of the mean component.
#' @param mu numerical vector containing the conditional mean of the response variable.
#' @param eta numerical vector containing the non linear predictor of the response variable.
#' @param phi numerical vector containing the conditional dispersion of the response variable.
#' @param f_mu function with signature `f(X, theta)` to compute the non linear predictor where `X` is the matrix containing the covariates and `theta` is the vector of parameters both of the mean component.
#' @param J_mu (optional) function to compute the jacobian of `f_mu()`. Should be a function with signature `J_mu(X, theta)` and should return a \eqn{n \times k} matrix where `n` is equal to `length(y)` and `k` is equal to the number of parameters of the mean component.
#' @param mu.eta function which computes the derivative of the mean component with respect to the non linear predictor as created by [make_link()] function.
#' @param variance function which computes the variance of the given family as a function of `mu` and `phi`.
#'
#' @return a numerical vector of length equal to the number of parameters for the mean component.
#' @export
#'
#' @examples
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
grad_mu <- function(y, X, beta, mu, eta, phi, f_mu, J_mu, mu.eta, variance){
  if(missing(J_mu)) J_mu <- makeJacobian(f_mu)
  w <- mu.eta(eta)/variance(mu, phi)
  j <- J_mu(X, beta)
  drop(t(w*j)%*%(y - mu))
}



#-------------------------------------------------------------------------------



#' Compute hessian with respect to parameters of mean component for Generalized Non-Linear Models
#'
#' @param y numerical vector of response variable.
#' @param X matrix containing the covariates for the mean component.
#' @param beta parameters of the mean component.
#' @param mu numerical vector containing the conditional mean of the response variable.
#' @param eta numerical vector containing the non linear predictor of the response variable.
#' @param phi numerical vector containing the conditional dispersion of the response variable.
#' @param f_mu function with signature `f(X, theta)` to compute the non linear predictor where `X` is the matrix containing the covariates and `theta` is the vector of parameters both of the mean component.
#' @param J_mu (optional) function to compute the jacobian of `f_mu()`. Should be a function with signature `J_mu(X, theta)` and should return a \eqn{n \times k} matrix where `n` is equal to `length(y)` and `k` is equal to the number of parameters of the mean component.
#' @param H_mu (optional) function to compute the hessian of `f_mu()`. Should be a function with signature `H_mu(X, theta)` and should return a list of \eqn{n} elements where `n` is equal to `length(y)` and each element should be the hessian of `f_mu()` computed at every row of `X`.
#' @param mu.eta function which computes the derivative of the mean component with respect to the non linear predictor as created by [make_link()] function.
#' @param mu2.eta2 function which computes the second derivative of the mean component with respect to the non linear predictor as created by [make_link()] function.
#' @param variance function which computes the variance of the given family as a function of `mu` and `phi`.
#'
#' @return a square symmetric matrix of dimension equal to the number of parameters for the mean component.
#' @export
#'
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
#'         fam_gamma$mu.eta, fam_gamma$mu2.eta2, fam_gamma$variance)
hess_mu <- function(y, X, beta, mu, eta, phi, f_mu, J_mu, H_mu, mu.eta, mu2.eta2, variance){
  if(missing(J_mu)) J_mu <- makeJacobian(f_mu)
  if(missing(H_mu)) H_mu <- makeHessian(f_mu)

  mu_eta <- mu.eta(eta)

  j <- J_mu(X, beta)
  v <- variance(mu, phi)
  w1 <- (mu_eta^2)/v
  h1 <- -crossprod(j, w1*j)


  w2 <- ((y - mu)/v)
  h21 <- crossprod(j, w2*mu2.eta2(eta)*j)

  h22 <- Reduce(`+`, Map(function(w, h) w*h, w = w2*mu_eta, h = H_mu(X, beta)))
  h1 + h21 + h22

}



#-------------------------------------------------------------------------------



#' Build full hessian of a Generalized Non-Linear model
#'
#' @param h_mu hessian matrix for the mean component.
#' @param h_phi hessian matrix for the dispersion component.
#' @param h_mu_phi hessian matrix for the mean and dispersion component.
#'
#' @return a \eqn{(k + p) \times (k + p)} square symmetric matrix where \eqn{k} and \eqn{p} are the number of parameters for the mean component and the dispersion component respectively.
#' @export
#'
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
  # if(expected){
  #   h_mu_phi <- matrix(0, nrow(h_mu_phi), ncol(h_mu_phi))
  # }
  rbind(cbind(h_mu, h_mu_phi), cbind(t(h_mu_phi), h_phi))
}



#-------------------------------------------------------------------------------



#' Compute minus log-likelihood for given values of parameters and a given log-likelihood function of Generalized Non-Linear Models
#'
#' @param y numerical vector of response variable.
#' @param X matrix containing the covariates for the mean component.
#' @param Z matrix containing the covariates for the dispersion component.
#' @param par numerical vector containing the parameters to evaluate the `loglik()` function at. It contains the parameters for both the mean and dispersion component.
#' @param npar_mu integer indicating the number of parameters of the mean component. It is used in order to differentiate between parameters for mean and dispersion component in `par`.
#' @param loglik a loglikelihood function to evaluate as created my a `gnlm_family()` function.
#' @param linkinv_mu inverse-link function for the mean component as created by [make_link()].
#' @param linkinv_phi inverse-link function for the dispersion component as created by [make_link()].
#'
#' @return numerical value of minus log-likelihood evaluated at given values of `par`.
#' @export
#'
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
#' log_loss(y, X, Z, par, npar_mu, fam_gamma$loglik, fam_gamma$linkinv_mu, fam_gamma$linkinv_phi)
log_loss <- function(y, X, Z, par, npar_mu, loglik, linkinv_mu, linkinv_phi){
  beta <- par[1:npar_mu]
  gamma <- par[(npar_mu + 1):length(par)]
  mu <- linkinv_mu(f_mu(X, beta))
  phi <- linkinv_phi(f_phi(Z, gamma))
  -sum(loglik(y, mu, phi))
}



#-------------------------------------------------------------------------------



