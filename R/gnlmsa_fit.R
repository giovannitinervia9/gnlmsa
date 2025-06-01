


#-------------------------------------------------------------------------------



#' Control Parameters for Newton-Raphson Algorithm
#'
#' Defines a set of user-configurable options to control the behavior of the Newton-Raphson
#' routine used in `gnlmsa` for optimization of the log-likelihood Generalized Non-Linear Models.
#'
#'
#' @param maxit Integer. Maximum number of iterations (default: 1000).
#' @param tol Numeric. Convergence tolerance on parameter updates (default: 1e-5).
#' @param expected Logical; if \code{TRUE}, the expected Fisher information is used in the optimization process.
#'   If \code{FALSE} (default), the observed Hessian is used.
#' @param unconstrained Logical; if \code{TRUE}, the optimization is performed in the unconstrained space (for detail see [`make_map_function()`] and [`reparametrize()`].
#'   If \code{FALSE} (default), the observed Hessian is the original parameter space.
#' @param verbose Logical; if \code{TRUE} (default), progress and warnings are printed.
#' @param regularization Numeric. A (small) constant added to the diagonal of the Hessian if it results non invertible (default: 1e-03).
#' @return A list containing the control parameters to be used in the Newton-Raphson routine.
#'
#'
#' @export
nr_control <- function(maxit = 100,
                       tol = 1e-05,
                       expected = FALSE,
                       unconstrained = TRUE,
                       verbose = TRUE,
                       regularization = 1e-03) {

  list(
    maxit = maxit,
    tol = tol,
    expected = expected,
    unconstrained = unconstrained,
    verbose = verbose,
    regularization = regularization
  )

}



#-------------------------------------------------------------------------------



#' Fit Generalized Non-Linear Models via Newton–Raphson Optimization
#'
#' Estimates the parameters of a Generalized Non-Linear Model (GNLM) by a full
#' Newton–Raphson (NR) routine.  The algorithm jointly updates the mean and
#' dispersion components using the gradient (score) and Hessian of the
#' log-likelihood. Arbitrary non-linear predictors may be supplied for both
#' components, with optional user-provided Jacobians and Hessians to improve
#' convergence.
#'
#' @section Details:
#' Iterations stop when the absolute change in every parameter is below
#' `nr_control_params$tol` or when `nr_control_params$maxit` iterations have
#' been reached.  If the log-likelihood decreases between successive iterations
#' a warning is issued, but the algorithm continues.  All distribution-specific
#' functions (link, variance, score, Hessians, log-likelihood) must be provided
#' through a [`family_gnlmsa`] object.
#'
#' The behaviour of the NR algorithm (maximum iterations, convergence
#' tolerance, whether to work in an unconstrained space, etc.) is governed by
#' the list returned by [`nr_control()`].  See `?nr_control` for a full
#' description of each option.
#'
#' @param y Numeric vector of response values.
#' @param X Design matrix for the mean component.
#' @param Z Design matrix for the dispersion component.
#' @param family A [`family_gnlmsa`] object defining the distribution,
#' link functions and associated derivatives.
#' @param f_mu Predictor function for the mean component.
#' @param J_mu (Optional) Jacobian of \code{f_mu}. Should return an \eqn{n \times k} matrix.
#' @param H_mu (Optional) Hessian of \code{f_mu}. Should return a list of \eqn{n} Hessian matrices.
#' @param f_phi Predictor function for the dispersion component.
#' @param J_phi (Optional) Jacobian of \code{f_phi}.
#' @param H_phi (Optional) Hessian of \code{f_phi}.
#' @param beta_start Initial values for parameters of the mean component.
#' @param lower_mu Numeric vector defining lower bounds for the mean parameters.
#' @param upper_mu Numeric vector defining upper bounds for the mean parameters.
#' @param gamma_start Initial values for parameters of the dispersion component.
#' @param lower_phi Numeric vector defining lower bounds for the dispersion parameters.
#' @param upper_phi Numeric vector defining upper bounds for the dispersion parameters.
#' @param fixed_params Optional list of two numeric vectors.  The first gives
#' the indices of parameters to be fixed, the second the
#' corresponding fixed values.
#' @param nr_control_params
#'                     A list of control settings created by
#'                     [`nr_control()`].  Defaults to `nr_control()`.  The list
#'                     elements are:
#'                     \itemize{
#'                       \item `maxit` – maximum NR iterations.
#'                       \item `tol` – convergence tolerance.
#'                       \item `expected` – logical; use expected (Fisher)
#'                             information when `TRUE`.
#'                       \item `unconstrained` – logical; optimise in the
#'                             unconstrained space when `TRUE`.
#'                       \item `verbose` – logical; print progress messages.
#'                       \item `regularization` – diagonal ridge added when the
#'                             Hessian is singular.
#'                     }
#'
#' @return An object of class `"gnlmsa_fit"` with components:
#' \describe{
#'   \item{beta}{Estimated mean-component coefficients.}
#'   \item{gamma}{Estimated dispersion-component coefficients.}
#'   \item{loglik}{Final log-likelihood.}
#'   \item{eta, mu}{Mean predictor and fitted values.}
#'   \item{vi, phi}{Dispersion predictor and fitted values.}
#'   \item{it}{Number of iterations used.}
#'   \item{maxit, tol}{Copied from `nr_control_params`.}
#'   \item{gradient, hessian}{Gradient and Hessian at the estimated values of the parameters (free
#'         parameters only).}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' y <- productivity$Y
#' X <- as.matrix(productivity[, -1])
#' Z <- cbind(1, prcomp(scale(X))$x[, 1])
#'
#' # Simulated–annealing warm start
#' sa <- sa_fit(
#'   y = y, X = X, Z = Z,
#'   family   = gnlmsa_Gamma("identity"),
#'   f_mu     = cobb_douglas()$f, J_mu = cobb_douglas()$J, H_mu = cobb_douglas()$H,
#'   f_phi    = Linear()$f,      J_phi = Linear()$J,      H_phi = Linear()$H,
#'   beta_start  = rep(1, 3), gamma_start = c(10, 0),
#'   mult = 1, nsim = 1,
#'   sa_control = sa_control(10000),
#'   expected   = TRUE, verbose = TRUE
#' )
#'
#' gnlmsa_fit(
#'   y = y, X = X, Z = Z,
#'   family = gnlmsa_Gamma("identity"),
#'   f_mu   = cobb_douglas()$f, J_mu = cobb_douglas()$J, H_mu = cobb_douglas()$H,
#'   f_phi  = Linear()$f,       J_phi = Linear()$J,      H_phi = Linear()$H,
#'   beta_start  = sa$beta, gamma_start = sa$gamma,
#'   nr_control_params = nr_control(maxit = 100, tol = 1e-5, expected = TRUE)
#' )
#' }
gnlmsa_fit <- function(y, X, Z, family,
                       f_mu, J_mu, H_mu,
                       f_phi, J_phi, H_phi,
                       beta_start, lower_mu, upper_mu,
                       gamma_start, lower_phi, upper_phi,
                       fixed_params = NULL,
                       nr_control_params = nr_control()) {

  if(missing(J_mu)) J_mu <- make_jacobian(f_mu)
  if(missing(H_mu)) H_mu <- make_hessian(f_mu)
  if(missing(J_phi)) J_phi <- make_jacobian(f_phi)
  if(missing(H_phi)) H_phi <- make_hessian(f_phi)

  maxit <- nr_control_params$maxit
  tol <- nr_control_params$tol
  expected <- nr_control_params$expected
  unconstrained <- nr_control_params$unconstrained
  verbose <- nr_control_params$verbose
  reg <- nr_control_params$regularization

  one_parameter <- family$family %in% c("gnlmsa_poisson", "gnlmsa_binomial")

  linkfun_mu <- family$linkfun_mu
  linkinv_mu <- family$linkinv_mu
  mu.eta <- family$mu.eta
  mu2.eta2 <- family$mu2.eta2

  linkfun_phi <- family$linkfun_phi
  linkinv_phi <- family$linkinv_phi
  phi.vi <- family$phi.vi
  phi2.vi2 <- family$phi2.vi2

  variance <- family$variance

  loglik <- family$loglik
  grad_mu <- family$grad_mu
  hess_mu <- family$hess_mu
  grad_phi <- family$grad_phi
  hess_phi <- family$hess_phi
  hess_mu_phi <- family$hess_mu_phi

  beta0 <- beta_start
  gamma0 <- gamma_start
  par0 <- c(beta0, gamma0)
  npar_mu <- length(beta0)
  npar <- length(par0)

  if (!is.null(fixed_params)) {
    fixed <- TRUE
    fixed_position <- fixed_params[[1]]
    fixed_values <- fixed_params[[2]]

    if (one_parameter) {
      fixed_position <- c(fixed_position, npar_mu + 1)
      fixed_values <- c(fixed_values, 1)
    }

  } else {

    if (one_parameter) {
      fixed <- TRUE
      fixed_position <- npar_mu + 1
      fixed_values <- 1
    } else {
      fixed <- FALSE
      fixed_position <- integer(0)
      fixed_values <- numeric(0)
    }
  }

  free_par <- setdiff(seq_len(npar), fixed_position)
  npar_free <- length(free_par)

  par0[fixed_position] <- fixed_values
  beta0 <- par0[1:npar_mu]
  gamma0 <- par0[(npar_mu + 1):npar]

  eta0 <- drop(f_mu(X, beta0))
  mu0 <- linkinv_mu(eta0)

  vi0 <- drop(f_phi(Z, gamma0))
  phi0 <- linkinv_phi(vi0)

  l0 <- sum(loglik(y, mu0, phi0))
  dev <- rep(2, length(par0))
  it <- 0
  par1 <- numeric(npar)
  par1[fixed_position] <- fixed_values

  lower <- c(lower_mu, lower_phi)[free_par]
  upper <- c(upper_mu, upper_phi)[free_par]
  map_functions <- make_map_function(lower, upper)
  map <- map_functions$map
  invert <- map_functions$invert
  jacobian <- map_functions$map_jacobian
  l_hist <- l0

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

  if (verbose) {
    cat("\n")
    cat("Newton-Raphson optimization:\n")
    cat("iter: ", 0,
        ", loglik: ", l0,
        ", mean abs. grad.: ", mean(abs(g[free_par])),
        "\n",
        sep = "")
  }

  while ((any(dev > tol)) && it <= maxit) {

    it <- it + 1

    if (unconstrained){
      rr <- reparametrize(par0, g, h, map_functions)
      step <- tryCatch(solve(rr$h_map[free_par, free_par], rr$g_map[free_par]), error = function(e) {
        solve(rr$h_map[free_par, free_par] + diag(reg, npar_free), rr$g_map[free_par])})
      par1[free_par] <- invert(map(par0[free_par]) - step)

    } else {
      step <- tryCatch(solve(h[free_par, free_par], g[free_par]), error = function(e) {
        solve(h[free_par, free_par] + diag(reg, npar_free), g[free_par])})
      par1[free_par] <- par0[free_par] - step
    }

    beta1 <- par1[1:npar_mu]
    gamma1 <- par1[(npar_mu + 1):npar]

    eta1 <- drop(f_mu(X, beta1))
    mu1 <- linkinv_mu(eta1)

    vi1 <- drop(f_phi(Z, gamma1))
    phi1 <- linkinv_phi(vi1)

    l1 <- sum(loglik(y, mu1, phi1))

    dev <- abs(par1 - par0)
    check <- l1 - l0 > 0

    par0 <- par1
    beta0 <- beta1
    gamma0 <- gamma1
    eta0 <- eta1
    mu0 <- mu1
    vi0 <- vi1
    phi0 <- phi1
    l0 <- l1
    l_hist <- c(l_hist, l1)

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

    if (verbose) {
      cat("iter: ", it,
          ", loglik: ", l0, if(!check) " (decreased!)",
          ", mean abs. grad.: ", mean(abs(g[free_par])),
          "\n",
          sep = "")
    }

  }

  out <- list(beta = beta0, gamma = gamma0, loglik = l0,
              eta = eta0, mu = mu0, vi = vi0, phi = phi0,
              it = it, maxit = maxit, tol = tol, l_hist = l_hist,
              gradient = g,
              hessian = h)
  class(out) <- "gnlmsa_fit"
  out

}

