


#-------------------------------------------------------------------------------




#' Family Objects for Generalized Non-Linear Models with Separate Mean and Dispersion Components
#'
#' Objects of class `family_gnlmsa` provide the necessary components to specify the distributional assumptions
#' and link functions for Generalized Non-Linear Models fitted via the [gnlmsa()] function.
#'
#' Unlike standard `family` objects used in generalized linear models (GLMs), `family_gnlmsa` accommodates
#' models where both the conditional mean \eqn{\mu} and dispersion \eqn{\phi} depend on covariates via
#' potentially nonlinear and distinct link functions.
#'
#' Additionally, the variance function is allowed to depend on both \eqn{\mu} and \eqn{\phi}, and
#' the log-likelihood gradient and Hessian with respect to dispersion parameters must be specified separately
#' for each distribution, as they are not generally reusable across families.
#'
#' @param object A `family_gnlmsa` object.
#' @param ... Additional arguments passed to methods.
#'
#' @returns A list of class `family_gnlmsa` with the following components:
#' \describe{
#'   \item{family}{Name of the assumed response distribution.}
#'
#'   \item{link_mu}{Name of the link function used for the mean component \eqn{\mu}.}
#'   \item{linkfun_mu}{The link function mapping \eqn{\mu} to the predictor \eqn{\eta}.}
#'   \item{linkinv_mu}{The inverse link function mapping \eqn{\eta} to \eqn{\mu}.}
#'   \item{mu.eta}{Derivative of \eqn{\mu} with respect to \eqn{\eta}.}
#'   \item{mu2.eta2}{Second derivative of \eqn{\mu} with respect to \eqn{\eta}.}
#'   \item{validmu}{Function to check that \eqn{\mu} values are within the valid parameter space.}
#'
#'   \item{link_phi}{Name of the link function used for the dispersion component \eqn{\phi}.}
#'   \item{linkfun_phi}{The link function mapping \eqn{\phi} to the predictor \eqn{v}.}
#'   \item{linkinv_phi}{The inverse link function mapping \eqn{v} to \eqn{\phi}.}
#'   \item{phi.vi}{Derivative of \eqn{\phi} with respect to \eqn{v}.}
#'   \item{phi2.vi2}{Second derivative of \eqn{\phi} with respect to \eqn{v}.}
#'   \item{validphi}{Function to check that \eqn{\phi} values are within the valid parameter space.}
#'
#'   \item{variance}{Function defining the conditional variance as a function of \eqn{\mu} and \eqn{\phi}.}
#'   \item{loglik}{Function that computes the log-likelihood contribution of each observation.}
#'   \item{grad_mu}{Gradient of the log-likelihood with respect to the mean parameters.}
#'   \item{hess_mu}{Hessian of the log-likelihood with respect to the mean parameters.}
#'   \item{grad_phi}{Gradient of the log-likelihood with respect to the dispersion parameters.}
#'   \item{hess_phi}{Hessian of the log-likelihood with respect to the dispersion parameters.}
#'   \item{hess_mu_phi}{Cross-partial derivative of the log-likelihood with respect to mean and dispersion parameters.}
#' }
#'
#' @export
#'
#' @examples
#' gnlmsa_Gamma(link_mu = "inverse", link_phi = "log")
family_gnlmsa <- function (object, ...) {
  UseMethod("family_gnlmsa")
}




#-------------------------------------------------------------------------------



#' Print Method for `family_gnlmsa` Objects
#'
#' Displays a summary of the distribution family and associated link functions used for modeling the mean and dispersion components in a `gnlmsa` model.
#'
#' This method is automatically invoked when printing a `family_gnlmsa` object.
#'
#' @param x An object of class `family_gnlmsa`.
#' @param ... Additional arguments passed to or from other methods (currently unused).
#'
#' @return The input object is returned invisibly. The function is called for its side effect of printing information to the console.
#' @export
#'
#' @examples
#' gamma_fam <- gnlmsa_Gamma()
#' gamma_fam
print.family_gnlmsa <- function(x, ...){
  cat("\nFamily:", x$family, "\n")
  cat("Link function:", x$link_mu, "\n")
  cat("Link function for dispersion:", x$link_phi, "\n\n")
  invisible(x)
}



#-------------------------------------------------------------------------------



#' Gamma Family for Generalized Non-Linear Models
#'
#' Constructs a `family_gnlmsa` object for use with Gamma-distributed response variables in
#' Generalized Non-Linear Models (GNLMs). This family supports separate link functions for both
#' the conditional mean \eqn{\mu} and the dispersion parameter \eqn{\phi}.
#'
#' @param link_mu A character string specifying the link function for the mean component \eqn{\mu}.
#'   Must be one of `"inverse"` (default), `"log"`, `"identity"` or `"sqrt"`.
#' @param link_phi A character string specifying the link function for the dispersion component \eqn{\phi}.
#'   Must be one of `"log"` (default), `"sqrt"` or `"identity"`.
#'
#' @return An object of class `family_gnlmsa` (inheriting from `family`) with the following components:
#' \describe{
#'   \item{family}{Name of the family: `"gnlmsa_Gamma"`}
#'
#'   \item{link_mu}{Name of the link function used for the mean.}
#'   \item{linkfun_mu}{Function mapping \eqn{\mu} to the predictor \eqn{\eta}.}
#'   \item{linkinv_mu}{Inverse of the link function, mapping \eqn{\eta} to \eqn{\mu}.}
#'   \item{mu.eta}{First derivative of \eqn{\mu} with respect to \eqn{\eta}.}
#'   \item{mu2.eta2}{Second derivative of \eqn{\mu} with respect to \eqn{\eta}.}
#'   \item{validmu}{Validator function for the mean: checks positivity and finiteness.}
#'
#'   \item{link_phi}{Name of the link function used for the dispersion.}
#'   \item{linkfun_phi}{Function mapping \eqn{\phi} to the predictor \eqn{v}.}
#'   \item{linkinv_phi}{Inverse of the link function, mapping \eqn{v} to \eqn{\phi}.}
#'   \item{phi.vi}{First derivative of \eqn{\phi} with respect to \eqn{v}.}
#'   \item{phi2.vi2}{Second derivative of \eqn{\phi} with respect to \eqn{v}.}
#'   \item{validphi}{Validator function for the dispersion: checks positivity and finiteness.}
#'
#'   \item{variance}{Function defining the conditional variance: \eqn{Var(Y_i \mid \mu_i, \phi_i) = \mu_i^2 / \phi_i}.}
#'   \item{loglik}{Log-likelihood contribution function for each observation.}
#'   \item{grad_mu}{Gradient of the log-likelihood with respect to the mean parameters.}
#'   \item{hess_mu}{Hessian of the log-likelihood with respect to the mean parameters.}
#'   \item{grad_phi}{Gradient of the log-likelihood with respect to the dispersion parameters.}
#'   \item{hess_phi}{Hessian of the log-likelihood with respect to the dispersion parameters.}
#'   \item{hess_mu_phi}{Cross-derivative of the log-likelihood with respect to mean and dispersion parameters.}
#'
#'   \item{simulate}{Function to simulate to be called by method `simulate`.}
#' }
#'
#' @details
#' The Gamma distribution is commonly used to model strictly positive, continuous, and skewed data.
#' This implementation assumes the following parameterization of the Gamma density:
#' \deqn{
#' f(y_i; \mu_i, \phi_i) =
#' \dfrac{ \left( \frac{\phi_i}{\mu_i} \right)^{\phi_i} y_i^{\phi_i - 1} \exp\left(- \phi_i \frac{y_i}{\mu_i} \right) }
#' {\Gamma(\phi_i)}
#' }
#'
#' The log-likelihood function (summed over all observations) is:
#' \deqn{
#' \ell(\mu, \phi) = \sum_{i=1}^n
#' \left\{
#' \phi_i \log \phi_i - \log \Gamma(\phi_i) + \phi_i \left( \log y_i - \log \mu_i - \frac{y_i}{\mu_i} \right) - \log y_i
#' \right\}
#' }
#'
#' The model allows for the mean and dispersion to depend on separate covariates via nonlinear predictors:
#' \deqn{\mu_i = g^{-1}( \eta(x_i, \beta) ), \quad \phi_i = h^{-1}( v(z_i, \gamma) )}
#'
#' where \eqn{g^{-1}} and \eqn{h^{-1}} are the inverse link functions for the mean and dispersion, respectively.
#'
#' @examples
#' # Gamma family with default links
#' gamma_fam <- gnlmsa_Gamma()
#'
#' # Gamma family with custom links
#' gamma_fam2 <- gnlmsa_Gamma(link_mu = "log", link_phi = "sqrt")
#'
#' @seealso [family_gnlmsa()], [gnlmsa()].
#' @importFrom stats rgamma
#' @export
gnlmsa_Gamma <- function(link_mu = "inverse", link_phi = "log"){

  linktemp_mu <- substitute(link_mu)
  linktemp_phi <- substitute(link_phi)

  okLinks_mu <- c("inverse", "log", "sqrt", "identity")
  okLinks_phi <- c("log", "sqrt", "identity")

  family <- "gnlmsa_Gamma"

  # create link for mu
  if (!is.character(linktemp_mu)) {linktemp_mu <- deparse(linktemp_mu)}
  if (linktemp_mu %in% okLinks_mu) {
    stats_mu <- make_link(linktemp_mu)
  }
  else if (is.character(linktemp_mu)) {
    stats_mu <- make_link(linktemp_mu)
    linktemp_mu <- link_mu
  }
  else {
    if (inherits(link_mu, "link-glm")) {
      stats_mu <- link_mu
      if (!is.null(stats_mu$name))
        linktemp_mu <- stats_mu$name
    }
    else {
      stop(gettextf("link \"%s\" for mu is not available for %s family; available links are %s",
                    linktemp_mu, family, paste(sQuote(okLinks_mu), collapse = ", ")),
           domain = NA)
    }
  }


  # create link for phi
  if (!is.character(linktemp_phi)) {
    linktemp_phi <- deparse(linktemp_phi)
  }
  if (linktemp_phi %in% okLinks_phi) {
    stats_phi <- make_link(linktemp_phi)
  }
  else if (is.character(linktemp_phi)) {
    stats_phi <- make_link(linktemp_phi)
    linktemp_phi <- link_phi
  }
  else {
    if (inherits(link_phi, "link-glm")) {
      stats_phi <- link_phi
      if (!is.null(stats_phi$name))
        linktemp_phi <- stats_phi$name
    }
    else {
      stop(gettextf("link \"%s\" for phi is not available for %s family; available links are %s",
                    linktemp_phi, family, paste(sQuote(okLinks_phi), collapse = ", ")),
           domain = NA)
    }
  }


  linkfun_mu <- stats_mu$linkfun
  linkinv_mu <- stats_mu$linkinv
  mu.eta <- stats_mu$mu.eta
  mu2.eta2 <- stats_mu$mu2.eta2

  linkfun_phi <- stats_phi$linkfun
  linkinv_phi <- stats_phi$linkinv
  phi.vi <- stats_phi$mu.eta
  phi2.vi2 <- stats_phi$mu2.eta2

  variance <- function(mu, phi) mu^2/phi
  validmu <- function(mu) all(is.finite(mu)) && all(mu > 0)
  validphi <- function(phi) all(is.finite(phi)) && all(phi > 0)
  loglik <- function(y, mu, phi){
    phi*log(phi) - lgamma(phi) + phi*(log(y) - log(mu)- y/mu) - log(y)
  }

  grad_mu <- ef_grad_mu

  hess_mu <- ef_hess_mu

  grad_phi <- function(y, Z, gamma, phi, vi, mu, f_phi, J_phi, phi.vi){
    if(missing(J_phi)) J_phi <- make_jacobian(f_phi)
    w <- (log(phi) + 1 - digamma(phi) + log(y) - log(mu) - y/mu)*phi.vi(vi)
    j <- J_phi(Z, gamma)
    colSums(j*w)
  }


  hess_phi <- function(y, Z, gamma, phi, vi, mu, f_phi, J_phi, H_phi,
                       phi.vi, phi2.vi2, expected = TRUE) {
    if (missing(J_phi)) J_phi <- make_jacobian(f_phi)
    if (missing(H_phi)) H_phi <- make_hessian(f_phi)

    j <- J_phi(Z, gamma)
    phi_vi <- phi.vi(vi)

    w1 <- (1 / phi - trigamma(phi)) * (phi_vi^2)
    h1 <- crossprod(j, w1 * j)

    if (expected) {
      return(h1)
    }

    # observed part
    w2 <- (log(phi) + 1 - digamma(phi) + log(y) - log(mu) - y / mu)
    h21 <- crossprod(j, j * w2 * phi2.vi2(vi))

    H_list <- H_phi(Z, gamma)
    h22 <- matrix(0, nrow = ncol(j), ncol = ncol(j))
    for (i in seq_along(w2)) {
      h22 <- h22 + w2[i] * phi_vi[i] * H_list[[i]]
    }

    h1 + h21 + h22
  }

  hess_mu_phi <- function(y, X, Z, beta, gamma,
                          mu, eta, phi, vi,
                          f_mu, J_mu, f_phi, J_phi,
                          mu.eta, phi.vi, expected = TRUE) {

    if (expected) {
      matrix(0, nrow = length(beta), ncol = length(gamma))
    } else {
      if (missing(J_mu)) {
        J_mu <- make_jacobian(f_mu)
      }
      if (missing(J_phi)) {
        J_phi <- make_jacobian(f_phi)
      }
      w <- ((y - mu)/mu^2)*mu.eta(eta)*phi.vi(vi)
      jmu <- J_mu(X, beta)
      jphi <- J_phi(Z, gamma)
      crossprod(w*jmu, jphi)
    }


  }

  simfun <- function(object, nsim) {
    mu <- object$mu
    phi <- object$phi
    rgamma(nsim*length(mu), phi, phi/mu)
  }


  out <- list(family = family,

              link_mu = linktemp_mu,
              linkfun_mu = linkfun_mu,
              linkinv_mu = linkinv_mu,
              mu.eta = mu.eta,
              mu2.eta2 = mu2.eta2,
              validmu = validmu,

              link_phi = linktemp_phi,
              linkfun_phi = linkfun_phi,
              linkinv_phi = linkinv_phi,
              phi.vi = phi.vi,
              phi2.vi2 = phi2.vi2,
              validphi = validphi,

              variance = variance,
              loglik = loglik,
              grad_mu = grad_mu,
              hess_mu = hess_mu,
              grad_phi = grad_phi,
              hess_phi = hess_phi,
              hess_mu_phi = hess_mu_phi,
              simfun = simfun)
  class(out) <- c("family_gnlmsa", "family")
  out
}



#-------------------------------------------------------------------------------



#' Gaussian Family for Generalized Non-Linear Models
#'
#' Constructs a `family_gnlmsa` object for use with Gaussian-distributed response variables in
#' Generalized Non-Linear Models (GNLMs). This family supports separate link functions for both
#' the conditional mean \eqn{\mu} and the dispersion parameter \eqn{\phi}.
#'
#' @param link_mu A character string specifying the link function for the mean component \eqn{\mu}.
#'   Must be one of `"inverse"` (default), `"log"`, `"identity"` or `"sqrt"`.
#' @param link_phi A character string specifying the link function for the dispersion component \eqn{\phi}.
#'   Must be one of `"log"` (default), `"sqrt"` or `"identity"`.
#'
#' @return An object of class `family_gnlmsa` (inheriting from `family`) with the following components:
#' \describe{
#'   \item{family}{Name of the family: `"gnlmsa_gaussian"`}

#'   \item{link_mu}{Name of the link function used for the mean.}
#'   \item{linkfun_mu}{Function mapping \eqn{\mu} to the predictor \eqn{\eta}.}
#'   \item{linkinv_mu}{Inverse of the link function, mapping \eqn{\eta} to \eqn{\mu}.}
#'   \item{mu.eta}{First derivative of \eqn{\mu} with respect to \eqn{\eta}.}
#'   \item{mu2.eta2}{Second derivative of \eqn{\mu} with respect to \eqn{\eta}.}
#'   \item{validmu}{Validator function for the mean: checks positivity and finiteness.}
#'
#'   \item{link_phi}{Name of the link function used for the dispersion.}
#'   \item{linkfun_phi}{Function mapping \eqn{\phi} to the predictor \eqn{v}.}
#'   \item{linkinv_phi}{Inverse of the link function, mapping \eqn{v} to \eqn{\phi}.}
#'   \item{phi.vi}{First derivative of \eqn{\phi} with respect to \eqn{v}.}
#'   \item{phi2.vi2}{Second derivative of \eqn{\phi} with respect to \eqn{v}.}
#'   \item{validphi}{Validator function for the dispersion: checks positivity and finiteness.}
#'
#'   \item{variance}{Function defining the conditional variance: \eqn{Var(Y_i \mid \mu_i, \phi_i) = \phi_i}.}
#'   \item{loglik}{Log-likelihood contribution function for each observation.}
#'   \item{grad_mu}{Gradient of the log-likelihood with respect to the mean parameters.}
#'   \item{hess_mu}{Hessian of the log-likelihood with respect to the mean parameters.}
#'   \item{grad_phi}{Gradient of the log-likelihood with respect to the dispersion parameters.}
#'   \item{hess_phi}{Hessian of the log-likelihood with respect to the dispersion parameters.}
#'   \item{hess_mu_phi}{Cross-derivative of the log-likelihood with respect to mean and dispersion parameters.}
#'
#'   \item{simulate}{Function to simulate to be called by method `simulate`.}
#' }
#'
#'
#' @seealso [family_gnlmsa()], [gnlmsa()].
#' @importFrom stats rnorm
#' @export
gnlmsa_gaussian <- function(link_mu = "identity", link_phi = "log"){

  linktemp_mu <- substitute(link_mu)
  linktemp_phi <- substitute(link_phi)

  okLinks_mu <- c("inverse", "log", "sqrt", "identity")
  okLinks_phi <- c("log", "sqrt", "identity")

  family <- "gnlmsa_gaussian"

  # create link for mu
  if (!is.character(linktemp_mu)) {linktemp_mu <- deparse(linktemp_mu)}
  if (linktemp_mu %in% okLinks_mu) {
    stats_mu <- make_link(linktemp_mu)
  }
  else if (is.character(linktemp_mu)) {
    stats_mu <- make_link(linktemp_mu)
    linktemp_mu <- link_mu
  }
  else {
    if (inherits(link_mu, "link-glm")) {
      stats_mu <- link_mu
      if (!is.null(stats_mu$name))
        linktemp_mu <- stats_mu$name
    }
    else {
      stop(gettextf("link \"%s\" for mu is not available for %s family; available links are %s",
                    linktemp_mu, family, paste(sQuote(okLinks_mu), collapse = ", ")),
           domain = NA)
    }
  }


  # create link for phi
  if (!is.character(linktemp_phi)) {
    linktemp_phi <- deparse(linktemp_phi)
  }
  if (linktemp_phi %in% okLinks_phi) {
    stats_phi <- make_link(linktemp_phi)
  }
  else if (is.character(linktemp_phi)) {
    stats_phi <- make_link(linktemp_phi)
    linktemp_phi <- link_phi
  }
  else {
    if (inherits(link_phi, "link-glm")) {
      stats_phi <- link_phi
      if (!is.null(stats_phi$name))
        linktemp_phi <- stats_phi$name
    }
    else {
      stop(gettextf("link \"%s\" for phi is not available for %s family; available links are %s",
                    linktemp_phi, family, paste(sQuote(okLinks_phi), collapse = ", ")),
           domain = NA)
    }
  }


  linkfun_mu <- stats_mu$linkfun
  linkinv_mu <- stats_mu$linkinv
  mu.eta <- stats_mu$mu.eta
  mu2.eta2 <- stats_mu$mu2.eta2

  linkfun_phi <- stats_phi$linkfun
  linkinv_phi <- stats_phi$linkinv
  phi.vi <- stats_phi$mu.eta
  phi2.vi2 <- stats_phi$mu2.eta2

  variance <- function(mu, phi) phi
  validmu <- function(mu) all(is.finite(mu))
  validphi <- function(phi) all(is.finite(phi)) && all(phi > 0)
  loglik <- function(y, mu, phi){
    -0.5*(log(2*pi) + log(phi) + (y - mu)^2/phi)
  }

  grad_mu <- ef_grad_mu

  hess_mu <- ef_hess_mu

  grad_phi <- function(y, Z, gamma, phi, vi, mu, f_phi, J_phi, phi.vi){
    if(missing(J_phi)) J_phi <- make_jacobian(f_phi)
    w <- (1/phi - ((y - mu)^2)/phi^2)*phi.vi(vi)
    j <- J_phi(Z, gamma)
    -0.5*colSums(j*w)
  }


  hess_phi <- function(y, Z, gamma, phi, vi, mu, f_phi, J_phi, H_phi,
                       phi.vi, phi2.vi2, expected = TRUE) {
    if (missing(J_phi)) J_phi <- make_jacobian(f_phi)
    if (missing(H_phi)) H_phi <- make_hessian(f_phi)

    j <- J_phi(Z, gamma)
    phi_vi <- phi.vi(vi)



    if (expected) {
      w1 <- (-1/phi^2 + 2/phi^2)*phi_vi^2
      h1 <- crossprod(j, w1 * j)
      return(-0.5*h1)
    } else {
      w1 <- (-1/phi^2 + (2*(y - mu)^2)/phi^3)*phi_vi^2
      h1 <- crossprod(j, w1 * j)
      # observed part
      w2 <- -1/phi^2 + (2*(y - mu)^2)/phi^3
      h21 <- crossprod(j, j * w2 * phi2.vi2(vi))

      H_list <- H_phi(Z, gamma)
      h22 <- matrix(0, nrow = ncol(j), ncol = ncol(j))
      for (i in seq_along(w2)) {
        h22 <- h22 + w2[i] * phi_vi[i] * H_list[[i]]
      }

      -0.5*(h1 + h21 + h22)
    }


  }

  hess_mu_phi <- function(y, X, Z, beta, gamma,
                          mu, eta, phi, vi,
                          f_mu, J_mu, f_phi, J_phi,
                          mu.eta, phi.vi, expected = TRUE) {

    if (expected) {
      matrix(0, nrow = length(beta), ncol = length(gamma))
    } else {
      if (missing(J_mu)) {
        J_mu <- make_jacobian(f_mu)
      }
      if (missing(J_phi)) {
        J_phi <- make_jacobian(f_phi)
      }
      w <- ((y - mu)/phi^2)*mu.eta(eta)*phi.vi(vi)
      jmu <- J_mu(X, beta)
      jphi <- J_phi(Z, gamma)
      crossprod(w*jmu, jphi)
    }


  }

  simfun <- function(object, nsim) {
    mu <- object$mu
    phi <- object$phi
    rnorm(nsim*length(mu), mu, sqrt(phi))
  }

  out <- list(family = family,

              link_mu = linktemp_mu,
              linkfun_mu = linkfun_mu,
              linkinv_mu = linkinv_mu,
              mu.eta = mu.eta,
              mu2.eta2 = mu2.eta2,
              validmu = validmu,

              link_phi = linktemp_phi,
              linkfun_phi = linkfun_phi,
              linkinv_phi = linkinv_phi,
              phi.vi = phi.vi,
              phi2.vi2 = phi2.vi2,
              validphi = validphi,

              variance = variance,
              loglik = loglik,
              grad_mu = grad_mu,
              hess_mu = hess_mu,
              grad_phi = grad_phi,
              hess_phi = hess_phi,
              hess_mu_phi = hess_mu_phi,
              simulate = simfun)
  class(out) <- c("family_gnlmsa", "family")
  out
}



#-------------------------------------------------------------------------------



#' Poisson Family for Generalized Non-Linear Models
#'
#' Constructs a `family_gnlmsa` object for use with Poisson-distributed response variables in
#' Generalized Non-Linear Models (GNLMs). In the Poisson the \eqn{\phi} parameter is fixed to be 1.
#'
#' @param link_mu A character string specifying the link function for the mean component \eqn{\mu}.
#'   Must be one of `"log"` (default), `"identity"` or `"sqrt"`.
#'
#' @return An object of class `family_gnlmsa` (inheriting from `family`) with the following components:
#' \describe{
#'   \item{family}{Name of the family: `"gnlmsa_poisson"`}
#'
#'   \item{link_mu}{Name of the link function used for the mean.}
#'   \item{linkfun_mu}{Function mapping \eqn{\mu} to the predictor \eqn{\eta}.}
#'   \item{linkinv_mu}{Inverse of the link function, mapping \eqn{\eta} to \eqn{\mu}.}
#'   \item{mu.eta}{First derivative of \eqn{\mu} with respect to \eqn{\eta}.}
#'   \item{mu2.eta2}{Second derivative of \eqn{\mu} with respect to \eqn{\eta}.}
#'   \item{validmu}{Validator function for the mean: checks positivity and finiteness.}
#'
#'   \item{link_phi}{Name of the link function used for the dispersion.}
#'   \item{linkfun_phi}{Function mapping \eqn{\phi} to the predictor \eqn{v}.}
#'   \item{linkinv_phi}{Inverse of the link function, mapping \eqn{v} to \eqn{\phi}.}
#'   \item{phi.vi}{First derivative of \eqn{\phi} with respect to \eqn{v}.}
#'   \item{phi2.vi2}{Second derivative of \eqn{\phi} with respect to \eqn{v}.}
#'   \item{validphi}{Validator function for the dispersion: checks positivity and finiteness.}
#'
#'   \item{variance}{Function defining the conditional variance: \eqn{Var(Y_i \mid \mu_i, \phi_i) = \mu_i}}
#'   \item{loglik}{Log-likelihood contribution function for each observation.}
#'   \item{grad_mu}{Gradient of the log-likelihood with respect to the mean parameters.}
#'   \item{hess_mu}{Hessian of the log-likelihood with respect to the mean parameters.}
#'   \item{grad_phi}{Gradient of the log-likelihood with respect to the dispersion parameters.}
#'   \item{hess_phi}{Hessian of the log-likelihood with respect to the dispersion parameters.}
#'   \item{hess_mu_phi}{Cross-derivative of the log-likelihood with respect to mean and dispersion parameters.}
#'
#'   \item{simulate}{Function to simulate to be called by method `simulate`.}
#' }
#'
#' @details
#' The Poisson distribution is commonly used to model count data.
#' The pdf of the Poisson distribution is:
#' \deqn{
#' f(y_i; \mu_i) =
#' \dfrac{\mu_i^{y_i} e^{\mu_i}}{y_i!}
#' }
#'
#' The log-likelihood function (summed over all observations) is:
#' \deqn{
#' \ell(\mu) = \sum_{i=1}^n
#' y_i \log \mu_i - \mu_i - \log y_i!
#' }
#'
#' @seealso [family_gnlmsa()], [gnlmsa()].
#' @importFrom stats rpois
#' @export
gnlmsa_poisson <- function(link_mu = "log"){

  link_phi <- "identity"

  linktemp_mu <- substitute(link_mu)
  linktemp_phi <- substitute(link_phi)

  okLinks_mu <- c("inverse", "log", "sqrt", "identity")
  okLinks_phi <- c("identity")

  family <- "gnlmsa_poisson"

  # create link for mu
  if (!is.character(linktemp_mu)) {linktemp_mu <- deparse(linktemp_mu)}
  if (linktemp_mu %in% okLinks_mu) {
    stats_mu <- make_link(linktemp_mu)
  }
  else if (is.character(linktemp_mu)) {
    stats_mu <- make_link(linktemp_mu)
    linktemp_mu <- link_mu
  }
  else {
    if (inherits(link_mu, "link-glm")) {
      stats_mu <- link_mu
      if (!is.null(stats_mu$name))
        linktemp_mu <- stats_mu$name
    }
    else {
      stop(gettextf("link \"%s\" for mu is not available for %s family; available links are %s",
                    linktemp_mu, family, paste(sQuote(okLinks_mu), collapse = ", ")),
           domain = NA)
    }
  }


  # create link for phi
  if (!is.character(linktemp_phi)) {
    linktemp_phi <- deparse(linktemp_phi)
  }
  if (linktemp_phi %in% okLinks_phi) {
    stats_phi <- make_link(linktemp_phi)
  }
  else if (is.character(linktemp_phi)) {
    stats_phi <- make_link(linktemp_phi)
    linktemp_phi <- link_phi
  }
  else {
    if (inherits(link_phi, "link-glm")) {
      stats_phi <- link_phi
      if (!is.null(stats_phi$name))
        linktemp_phi <- stats_phi$name
    }
    else {
      stop(gettextf("link \"%s\" for phi is not available for %s family; available links are %s",
                    linktemp_phi, family, paste(sQuote(okLinks_phi), collapse = ", ")),
           domain = NA)
    }
  }


  linkfun_mu <- stats_mu$linkfun
  linkinv_mu <- stats_mu$linkinv
  mu.eta <- stats_mu$mu.eta
  mu2.eta2 <- stats_mu$mu2.eta2

  linkfun_phi <- stats_phi$linkfun
  linkinv_phi <- stats_phi$linkinv
  phi.vi <- stats_phi$mu.eta
  phi2.vi2 <- stats_phi$mu2.eta2

  variance <- function(mu, phi = 1) mu
  validmu <- function(mu) all(is.finite(mu)) && all(mu > 0)
  validphi <- function(phi) all(is.finite(phi)) && all(phi > 0)
  loglik <- function(y, mu, phi = 1){
    y*log(mu) - mu - lgamma(y + 1)
  }

  grad_mu <- ef_grad_mu

  hess_mu <- ef_hess_mu

  grad_phi <- function(y, Z, gamma, phi = 1, vi, mu, f_phi, J_phi, phi.vi){
    NA
  }


  hess_phi <- function(y, Z, gamma, phi, vi, mu, f_phi, J_phi, H_phi,
                       phi.vi, phi2.vi2, expected = TRUE) {
    matrix(NA)
  }

  hess_mu_phi <- function(y, X, Z, beta, gamma,
                          mu, eta, phi, vi,
                          f_mu, J_mu, f_phi, J_phi,
                          mu.eta, phi.vi, expected = TRUE) {
    matrix(NA, length(beta), 1)
  }

  simfun <- function(object, nsim) {
    mu <- object$mu
    rpois(nsim*length(mu), mu)
  }


  out <- list(family = family,

              link_mu = linktemp_mu,
              linkfun_mu = linkfun_mu,
              linkinv_mu = linkinv_mu,
              mu.eta = mu.eta,
              mu2.eta2 = mu2.eta2,
              validmu = validmu,

              link_phi = linktemp_phi,
              linkfun_phi = linkfun_phi,
              linkinv_phi = linkinv_phi,
              phi.vi = phi.vi,
              phi2.vi2 = phi2.vi2,
              validphi = validphi,

              variance = variance,
              loglik = loglik,
              grad_mu = grad_mu,
              hess_mu = hess_mu,
              grad_phi = grad_phi,
              hess_phi = hess_phi,
              hess_mu_phi = hess_mu_phi,
              simfun = simfun)
  class(out) <- c("family_poisson", "family")
  out
}



#-------------------------------------------------------------------------------



#' Binomial Family for Generalized Non-Linear Models
#'
#' Constructs a `family_gnlmsa` object for use with Binomial-distributed response variables in
#' Generalized Non-Linear Models (GNLMs). In the Binomial model the \eqn{\phi} parameter is fixed to be 1.
#'
#' @param link_mu A character string specifying the link function for the mean component \eqn{\mu}.
#'
#' @return An object of class `family_gnlmsa` (inheriting from `family`) with the following components:
#' \describe{
#'   \item{family}{Name of the family: `"gnlmsa_binomial"`}
#'
#'   \item{link_mu}{Name of the link function used for the mean.}
#'   \item{linkfun_mu}{Function mapping \eqn{\mu} to the predictor \eqn{\eta}.}
#'   \item{linkinv_mu}{Inverse of the link function, mapping \eqn{\eta} to \eqn{\mu}.}
#'   \item{mu.eta}{First derivative of \eqn{\mu} with respect to \eqn{\eta}.}
#'   \item{mu2.eta2}{Second derivative of \eqn{\mu} with respect to \eqn{\eta}.}
#'   \item{validmu}{Validator function for the mean: checks positivity and finiteness.}
#'
#'   \item{link_phi}{Name of the link function used for the dispersion.}
#'   \item{linkfun_phi}{Function mapping \eqn{\phi} to the predictor \eqn{v}.}
#'   \item{linkinv_phi}{Inverse of the link function, mapping \eqn{v} to \eqn{\phi}.}
#'   \item{phi.vi}{First derivative of \eqn{\phi} with respect to \eqn{v}.}
#'   \item{phi2.vi2}{Second derivative of \eqn{\phi} with respect to \eqn{v}.}
#'   \item{validphi}{Validator function for the dispersion: checks positivity and finiteness.}
#'
#'   \item{variance}{Function defining the conditional variance: \eqn{Var(Y_i \mid \mu_i, \phi_i) = \mu_i(1 - \mu_i)}}
#'   \item{loglik}{Log-likelihood contribution function for each observation.}
#'   \item{grad_mu}{Gradient of the log-likelihood with respect to the mean parameters.}
#'   \item{hess_mu}{Hessian of the log-likelihood with respect to the mean parameters.}
#'   \item{grad_phi}{Gradient of the log-likelihood with respect to the dispersion parameters.}
#'   \item{hess_phi}{Hessian of the log-likelihood with respect to the dispersion parameters.}
#'   \item{hess_mu_phi}{Cross-derivative of the log-likelihood with respect to mean and dispersion parameters.}
#'
#'   \item{simulate}{Function to simulate to be called by method `simulate`.}
#' }
#'
#' @details
#' The Binomial distribution is commonly used to model count data.
#' The pdf of the Binomial distribution is:
#' \deqn{
#' f(y_i; \mu_i) =
#' \dbinom{m_i}{y_i} \mu_i^{y_i} (1 - \mu_i)^{m_i - y_i}
#' }
#'
#' The log-likelihood function (summed over all observations) is:
#' \deqn{
#' \ell(\mu) = \sum_{i=1}^n
#' \log \dbinom{m_i}{y_i} + y_i \log \mu_i + (m_i - y_i) \log (1 - \mu_i)
#' }
#'
#' @seealso [family_gnlmsa()], [gnlmsa()].
#' @importFrom stats rpois
#' @export
gnlmsa_binomial <- function(link_mu = "logit"){

  link_phi <- "identity"

  linktemp_mu <- substitute(link_mu)
  linktemp_phi <- substitute(link_phi)

  okLinks_mu <- c("logit", "probit", "cloglog", "cauchit", "log")
  okLinks_phi <- c("identity")

  family <- "gnlmsa_binomial"

  # create link for mu
  if (!is.character(linktemp_mu)) {linktemp_mu <- deparse(linktemp_mu)}
  if (linktemp_mu %in% okLinks_mu) {
    stats_mu <- make_link(linktemp_mu)
  }
  else if (is.character(linktemp_mu)) {
    stats_mu <- make_link(linktemp_mu)
    linktemp_mu <- link_mu
  }
  else {
    if (inherits(link_mu, "link-glm")) {
      stats_mu <- link_mu
      if (!is.null(stats_mu$name))
        linktemp_mu <- stats_mu$name
    }
    else {
      stop(gettextf("link \"%s\" for mu is not available for %s family; available links are %s",
                    linktemp_mu, family, paste(sQuote(okLinks_mu), collapse = ", ")),
           domain = NA)
    }
  }


  # create link for phi
  if (!is.character(linktemp_phi)) {
    linktemp_phi <- deparse(linktemp_phi)
  }
  if (linktemp_phi %in% okLinks_phi) {
    stats_phi <- make_link(linktemp_phi)
  }
  else if (is.character(linktemp_phi)) {
    stats_phi <- make_link(linktemp_phi)
    linktemp_phi <- link_phi
  }
  else {
    if (inherits(link_phi, "link-glm")) {
      stats_phi <- link_phi
      if (!is.null(stats_phi$name))
        linktemp_phi <- stats_phi$name
    }
    else {
      stop(gettextf("link \"%s\" for phi is not available for %s family; available links are %s",
                    linktemp_phi, family, paste(sQuote(okLinks_phi), collapse = ", ")),
           domain = NA)
    }
  }


  linkfun_mu <- stats_mu$linkfun
  linkinv_mu <- stats_mu$linkinv
  mu.eta <- stats_mu$mu.eta
  mu2.eta2 <- stats_mu$mu2.eta2

  linkfun_phi <- stats_phi$linkfun
  linkinv_phi <- stats_phi$linkinv
  phi.vi <- stats_phi$mu.eta
  phi2.vi2 <- stats_phi$mu2.eta2

  variance <- function(mu, phi = 1) mu*(1 - mu)
  validmu <- function(mu) all(mu > 0 & mu < 1)
  validphi <- function(phi) all(is.finite(phi)) && all(phi > 0)
  loglik <- function(y, mu, phi = 1){
    m <- 1
    lgamma(m + 1) - lgamma(y + 1) - lgamma(m - y + 1) + y*log(mu) + (m - y)*log(1 - mu)
  }

  grad_mu <- ef_grad_mu

  hess_mu <- ef_hess_mu

  grad_phi <- function(y, Z, gamma, phi = 1, vi, mu, f_phi, J_phi, phi.vi){
    NA
  }


  hess_phi <- function(y, Z, gamma, phi, vi, mu, f_phi, J_phi, H_phi,
                       phi.vi, phi2.vi2, expected = TRUE) {
    matrix(NA)
  }

  hess_mu_phi <- function(y, X, Z, beta, gamma,
                          mu, eta, phi, vi,
                          f_mu, J_mu, f_phi, J_phi,
                          mu.eta, phi.vi, expected = TRUE) {
    matrix(NA, length(beta), 1)
  }

  simfun <- function(object, nsim) {
    mu <- object$mu
    rpois(nsim*length(mu), mu)
  }


  out <- list(family = family,

              link_mu = linktemp_mu,
              linkfun_mu = linkfun_mu,
              linkinv_mu = linkinv_mu,
              mu.eta = mu.eta,
              mu2.eta2 = mu2.eta2,
              validmu = validmu,

              link_phi = linktemp_phi,
              linkfun_phi = linkfun_phi,
              linkinv_phi = linkinv_phi,
              phi.vi = phi.vi,
              phi2.vi2 = phi2.vi2,
              validphi = validphi,

              variance = variance,
              loglik = loglik,
              grad_mu = grad_mu,
              hess_mu = hess_mu,
              grad_phi = grad_phi,
              hess_phi = hess_phi,
              hess_mu_phi = hess_mu_phi,
              simfun = simfun)
  class(out) <- c("family_binomial", "family")
  out
}
