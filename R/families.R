


#-------------------------------------------------------------------------------



#' Print method for `family_gnlmsa` objects
#'
#' @param x A `family_gnlmsa` object
#' @param ... further arguments passed to methods.
#'
#' @return Details about the family definition and the family object.
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



#' Gamma family for Generalized Non-Linear Models
#'
#' This function creates a family object for Gamma distributions in Generalized Non-Linear Models (GNLM).
#'
#' @param link_mu Character string specifying the link function for the mean component.
#'   Must be one of "inverse" (default), "log", or "identity".
#' @param link_phi Character string specifying the link function for the dispersion component.
#'   Must be one of "log" (default) or "sqrt".
#'
#' @return A \code{family_gnlmsa} object (which also inherits from the \code{family} class) containing:
#' \describe{
#'   \item{family}{Name of the family distribution ("gnlmsa_Gamma")}
#'   \item{link_mu}{Name of the link function for \eqn{\mu}}
#'   \item{linkfun_mu}{Function that transforms \eqn{\mu} to \eqn{\eta}}
#'   \item{linkinv_mu}{Inverse link function that transforms \eqn{\eta} to \eqn{\mu}}
#'   \item{mu.eta}{First derivative \eqn{\partial \mu / \partial \eta}}
#'   \item{mu2.eta2}{Second derivative \eqn{\partial^2 \mu / \partial^2 \eta^2}}
#'   \item{validmu}{Function to check if mean values are valid}
#'
#'   \item{link_phi}{Name of the link function for \eqn{\phi}}
#'   \item{linkfun_phi}{Function that transforms \eqn{\phi} to \eqn{v}}
#'   \item{linkinv_phi}{Inverse link function that transforms \eqn{v} to \eqn{\phi}}
#'   \item{phi.vi}{First derivative \eqn{\partial \phi / \partial v}}
#'   \item{phi2.vi2}{Second derivative \eqn{\partial^2 \phi / \partial^2 v^2}}
#'   \item{validphi}{Function to check if dispersion values are valid}
#'
#'   \item{variance}{Function defining the mean-variance relationship: \eqn{\mu^2 / \phi}}
#'   \item{loglik}{Function to calculate the log-likelihood}
#'   \item{grad_phi}{Function for gradient of log-likelihood with respect to dispersion parameters}
#'   \item{hess_phi}{Function for Hessian matrix of log-likelihood with respect to dispersion parameters}
#'   \item{hess_mu_phi}{Function for cross-derivative Hessian between mean and dispersion parameters}
#' }
#'
#' @details
#' The Gamma distribution is commonly used for modeling continuous, positive data with right skew.
#' In this implementation, the pdf of the Gamma distribution is
#' \deqn{f(y_i, \phi_i, \mu_i) = \dfrac{ \left(\dfrac{\phi_i}{\mu_i}\right)^{\phi_i}  y_i^{\phi_i - 1} \exp \left(- \phi_i \dfrac{y_i}{\mu_i} \right)   }    {\Gamma(\phi_i)} }
#'
#' The log-likelihood function is
#' \deqn{\ell(\mu, \phi) = \sum_{i = 1}^n \left \{ \phi_i \log \phi_i - \log \Gamma(\phi_i) + \phi_i\left( \log y_i - \log \mu_i - \dfrac{y_i}{\mu_i} \right) - \log y_i \right\}}
#'
#' with
#' \deqn{\mu_i = g^{-1}(\eta (x_i, \beta))}
#' \deqn{\phi_i = h^{-1}(v(z_i, \gamma))}
#'
#' @examples
#' # Create a Gamma family with inverse link for mean and log link for dispersion
#' gamma_fam <- gnlmsa_Gamma()
#'
#' # Create a Gamma family with log link for mean and sqrt link for dispersion
#' gamma_fam2 <- gnlmsa_Gamma(link_mu = "log", link_phi = "sqrt")
#'
#' @seealso Other family functions in the gnlmsa package
#'
#' @export
gnlmsa_Gamma <- function(link_mu = "inverse", link_phi = "log"){

  linktemp_mu <- substitute(link_mu)
  linktemp_phi <- substitute(link_phi)

  okLinks_mu <- c("inverse", "log", "identity")
  okLinks_phi <- c("log", "sqrt")

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

  grad_phi <- function(y, Z, gamma, phi, vi, mu, f_phi, J_phi, phi.vi){
    if(missing(J_phi)) J_phi <- make_jacobian(f_phi)
    w <- (log(phi) + 1 - digamma(phi) + log(y) - log(mu) - y/mu)*phi.vi(vi)
    j <- J_phi(Z, gamma)
    colSums(j*w)
  }


  hess_phi <- function(y, Z, gamma, phi, vi, mu, f_phi, J_phi, H_phi, phi.vi, phi2.vi2, expected = TRUE){
    if(missing(J_phi)) J_phi <- make_jacobian(f_phi)
    if(missing(H_phi)) H_phi <- make_hessian(f_phi)

    phi_vi <- phi.vi(vi)

    j <- J_phi(Z, gamma)
    w1 <- (1/phi - trigamma(phi))*phi_vi^2
    h1 <- crossprod(j, w1*j)

    if (expected) {
      return(h1)
    } else {
      w2 <- (log(phi) + 1 - digamma(phi) + log(y) - log(mu) - y/mu)

      h21 <- crossprod(j, j*w2*phi2.vi2(vi))
      h22 <- Reduce(`+`, Map(function(w, h) w*h, w = w2*phi_vi, h = H_phi(Z, gamma)))
      h1 + h21 + h22
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
      w <- ((y - mu)/mu^2)*mu.eta(eta)*phi.vi(vi)
      jmu <- J_mu(X, beta)
      jphi <- J_phi(Z, gamma)
      crossprod(w*jmu, jphi)
    }


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
              grad_phi = grad_phi,
              hess_phi = hess_phi,
              hess_mu_phi = hess_mu_phi)
  class(out) <- c("family_gnlmsa", "family")
  out
}



#-------------------------------------------------------------------------------



