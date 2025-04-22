#' Fit Generalized Non‑Linear Models via Newton–Raphson Algorithm
#'
#' Simultaneously estimates the mean and dispersion components of a
#' Generalized Non‑Linear Model (GNLM) by a full Newton–Raphson routine that
#' relies on the score vector and observed (or, optionally, expected) Hessian
#' matrix.
#' Arbitrary non‑linear predictors can be supplied for both components, and
#' users may optionally provide analytic Jacobians and Hessians to accelerate
#' convergence.
#'
#' @section Details:
#' The iterative loop stops when the absolute change of every parameter falls
#' below `tol` or when `maxit` iterations are reached.
#' A warning is issued if the log‑likelihood decreases between successive
#' iterations, but the algorithm continues from the current values.
#' The required distributional functions (links, variance, score, Hessians and
#' log‑likelihood) must be provided through a [`family_gnlmsa`] object.
#'
#' @param y numerical vector of response values.
#' @param X design matrix for the mean component.
#' @param Z design matrix for the dispersion component.
#' @param family a [`family_gnlmsa`] object.
#' @param f_mu predictor function for the mean_component.
#' @param J_mu (optional) jacobian of `f_mu()`.
#' @param H_mu (optional) hessian of `f_mu()`.
#' @param f_phi predictor function for the dispersion component.
#' @param J_phi (optional) jacobian of `f_phi()`.
#' @param H_phi (optional) hessian of `f_phi()`.
#' @param beta_start numerical vector of starting values for the parameters of the mean component.
#' @param gamma_start numerical vector of starting values for the parameters of the dispersion component.
#' @param maxit integer. Maximum number of iterations.
#' @param tol convergence tolerance on the parameters (default 1e‑5).
#' @param expected logical indicating whether to use the expected (default, `expected = TRUE`) or the observed (`expected = FALSE`) in the optimization algorithm.
#'
#' @return
#' A list of class `"gnlmsa_fit"` with the following components:
#' \describe{
#'   \item{beta}{Estimated coefficients for the mean component.}
#'   \item{gamma}{Estimated coefficients for the dispersion component.}
#'   \item{loglik}{Final log‑likelihood value.}
#'   \item{eta}{Vector of linear predictors for the mean,
#'              \eqn{\hat{\eta} = \eta(X, \hat{\beta})}.}
#'   \item{mu}{Fitted means, \eqn{\hat{\mu} = g^{-1}(\hat{\eta})}.}
#'   \item{vi}{Vector of linear predictors for the dispersion,
#'              \eqn{\hat{v} = v(Z, \hat{\gamma})}.}
#'   \item{phi}{Fitted dispersion values, \eqn{\hat{\phi} = h^{-1}(\hat{v})}.}
#'   \item{it}{Number of iterations performed.}
#'   \item{maxit, tol}{The stopping criteria that were used.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' y <- productivity$Y
#' X <- as.matrix(productivity[, -1])
#' Z <- cbind(1, prcomp(scale(X))$x[, 1])
#'
#' sa <- sa_fit(y = y, X = X, Z = Z, family = gnlmsa_Gamma("identity"),
#'              f_mu = cobb_douglas()$f, J_mu = cobb_douglas()$J,
#'              H_mu = cobb_douglas()$H,
#'              f_phi = Linear()$f, J_phi = Linear()$J, H_phi = Linear()$H,
#'              beta_start = rep(1, 3), lower_mu = rep(0, 3), upper_mu = rep(Inf, 3),
#'              gamma_start = c(10, 0), lower_phi = rep(-Inf, 2), upper_phi = rep(Inf, 2),
#'              mult = 1, nsim = 1, sa_control = sa_control(10000), expected = TRUE,
#'              verbose = TRUE)
#'
#' gnlmsa_fit(y = y, X = X, Z = Z, family = gnlmsa_Gamma("identity"),
#'            f_mu = cobb_douglas()$f, J_mu = cobb_douglas()$J,
#'            H_mu = cobb_douglas()$H,
#'            f_phi = Linear()$f, J_phi = Linear()$J, H_phi = Linear()$H,
#'            beta_start = sa$beta, gamma_start = sa$gamma,
#'            maxit = 100, tol = 1e-05, expected = TRUE)
#' }
gnlmsa_fit <- function(y, X, Z, family,
                       f_mu, J_mu, H_mu,
                       f_phi, J_phi, H_phi,
                       beta_start, gamma_start, maxit = 100, tol = 1e-05,
                       expected = TRUE) {

  if(missing(J_mu)) J_mu <- make_jacobian(f_mu)
  if(missing(H_mu)) H_mu <- make_hessian(f_mu)
  if(missing(J_phi)) J_phi <- make_jacobian(f_phi)
  if(missing(H_phi)) H_phi <- make_hessian(f_phi)

  linkfun_mu <- family$linkfun_mu
  linkinv_mu <- family$linkinv_mu
  mu.eta <- family$mu.eta
  mu2.eta2 <- family$mu2.eta2

  linkfun_phi <- family$linkfun_phi
  linkinv_phi <- family$linkinv_phi
  phi.vi <- family$phi.vi
  phi2.vi2 <- family$phi2.vi2

  variance <- family$variance

  grad_phi <- family$grad_phi
  hess_phi <- family$hess_phi
  hess_mu_phi <- family$hess_mu_phi
  loglik <- family$loglik


  beta0 <- beta_start
  gamma0 <- gamma_start
  par0 <- c(beta0, gamma0)
  npar_mu <- length(beta0)

  eta0 <- f_mu(X, beta0)
  mu0 <- linkinv_mu(eta0)

  vi0 <- f_phi(Z, gamma0)
  phi0 <- linkinv_phi(vi0)

  l0 <- sum(loglik(y, mu0, phi0))
  dev <- rep(2, length(par0))
  it <- 0

  while ((any(dev > tol)) && it <= maxit) {

    it <- it + 1

    g_mu <- grad_mu(y, X, beta0, mu0, eta0, phi0, f_mu, J_mu, mu.eta, variance)
    g_phi <- grad_phi(y, Z, gamma0, phi0, vi0, mu0, f_phi, J_phi, phi.vi)
    g <- c(g_mu, g_phi)

    h_mu <- hess_mu(y, X, beta0, mu0, eta0, phi0, f_mu, J_mu, H_mu, mu.eta, mu2.eta2,
                    variance, expected)
    h_phi <- hess_phi(y, Z, gamma0, phi0, vi0, mu0, f_phi, J_phi, H_phi, phi.vi,
                      phi2.vi2, expected)
    h_mu_phi <- hess_mu_phi(y, X, Z, beta0, gamma0, mu0, eta0, phi0, vi0, f_mu, J_mu,
                            f_phi, J_phi, mu.eta, phi.vi, expected)
    h <- hess(h_mu, h_phi, h_mu_phi)

    par1 <- par0 - solve(h, g)
    beta1 <- par1[1:npar_mu]
    gamma1 <- par1[(npar_mu + 1):length(par1)]

    eta1 <- f_mu(X, beta1)
    mu1 <- linkinv_mu(eta1)

    vi1 <- f_phi(Z, gamma1)
    phi1 <- linkinv_phi(vi1)

    l1 <- sum(loglik(y, mu1, phi1))

    dev <- abs(par1 - par0)
    check <- l1 - l0 > 0
    if (!check) {
      warning("log-likelihood decreased at iteration ", it)
    }

    par0 <- par1
    beta0 <- beta1
    gamma0 <- gamma1
    eta0 <- eta1
    mu0 <- mu1
    vi0 <- vi1
    phi0 <- phi1
    l0 <- l1

  }

  out <- list(beta = beta0, gamma = gamma0, loglik = l0,
              eta = eta0, mu = mu0, vi = vi0, phi = phi0,
              it = it, maxit = maxit, tol = tol)
  class(out) <- "gnlmsa_fit"
  out

}

