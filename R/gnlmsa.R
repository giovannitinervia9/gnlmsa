#' Fitting Generalized Non-Linear Models
#'
#' @param y numerical vector of response values.
#' @param X design matrix for the mean model.
#' @param Z design matrix for the dispersion model.
#' @param family a [`family_gnlmsa`] object.
#' @param f_mu predictor function for the mean component.
#' @param J_mu (optional) Jacobian of `f_mu` with respect to the parameters.
#' @param H_mu (optional) Hessian of `f_mu` with respect to the parameters.
#' @param f_phi predictor function for the dispersion component.
#' @param J_phi (optional) Jacobian of `f_phi` with respect to the parameters.
#' @param H_phi (optional) Hessian of `f_phi` with respect to the parameters.
#' @param beta_start initial values for the parameters of the mean component.
#' @param lower_mu lower bounds for the parameters of the mean component.
#' @param upper_mu upper bounds for the parameters of the mean component.
#' @param gamma_start initial values for the parameters of the dispersion component.
#' @param lower_phi lower bounds for the parameters of the dispersion component.
#' @param upper_phi upper bounds for the parameters of the dispersion component.
#' @param mult numerical vector multipliers of errors in Simulated Annealing algorithm.
#' @param nsim number of Simulated Annealing simulations to perform (possibly in parallel).
#' @param sa_control list of control parameters for Simulated Annealing algorithm, created with `sa_control()`
#' @param maxit integer. Maximum number of iterations in Newton-Raphson optimization.
#' @param tol convergence tolerance on the parameters (default 1e‑5).
#' @param expected logical indicating whether to use the expected (default, `expected = TRUE`) or the observed (`expected = FALSE`) in the optimization algorithm.
#' @param verbose logical indicating whether to print progress information.
#' @param beta_names (optional) character vector containing the names of parameters of mean component. Must be of the same length of `beta_start`, `lower_mu` and `upper_mu`.
#' @param beta_names (optional) character vector containing the names of parameters of dispersion component. Must be of the same length of `gamma_start`, `lower_phi` and `upper_phi`.
#'
#' @returns A list.
#' @export
#'
#' @examples
#' \dontrun{
#' y <- productivity$Y
#' X <- as.matrix(productivity[, -1])
#' Z <- cbind(1, prcomp(scale(X))$x[, 1])
#' family <- gnlmsa_Gamma("identity")
#' f_mu <- cobb_douglas()$f
#' J_mu <- cobb_douglas()$J
#' H_mu <- cobb_douglas()$H
#' f_phi <- Linear()$f
#' J_phi <- Linear()$J
#' H_phi <- Linear()$H
#' beta_start <- c(10, .5, .5)
#' lower_mu <- rep(0, 3)
#' upper_mu <- rep(Inf, 3)
#' gamma_start <- c(log(var(y)), 0)
#' lower_phi <- rep(-Inf, 2)
#' upper_phi <- rep(Inf, 2)
#' mult <- 1
#' nsim <- 1
#' sa_control <- sa_control(10000)
#' maxit <- 100
#' tol <- 1e-05
#' expected <- TRUE
#' verbose <- TRUE
#' beta_names <- c("a", "alpha", "beta")
#' gnlmsa(y, X, Z, family,
#'        f_mu, J_mu, H_mu,
#'        f_phi, J_phi, H_phi,
#'        beta_start, lower_mu, upper_mu,
#'        gamma_start, lower_phi, upper_phi,
#'        mult, nsim, sa_control = sa_control(),
#'        maxit = 100, tol = 1e-05,
#'        expected = TRUE, verbose = TRUE,
#'        beta_names)
#' }
gnlmsa <- function (y, X, Z, family,
                    f_mu, J_mu, H_mu,
                    f_phi, J_phi, H_phi,
                    beta_start, lower_mu, upper_mu,
                    gamma_start, lower_phi, upper_phi,
                    mult, nsim, sa_control = sa_control(),
                    maxit = 100, tol = 1e-05,
                    expected = TRUE, verbose = TRUE,
                    beta_names, gamma_names) {


  if(length(lower_mu) != length(upper_mu)) stop("lower_mu and upper_mu must have the same length")
  if(length(lower_phi) != length(upper_phi)) stop("lower_phi and upper_phi must have the same length")

  if(missing(J_mu)) J_mu <- make_jacobian(f_mu)
  if(missing(H_mu)) H_mu <- make_hessian(f_mu)
  if(missing(J_phi)) J_phi <- make_jacobian(f_phi)
  if(missing(H_phi)) H_phi <- make_hessian(f_phi)

  X <- as.matrix(X)
  Z <- as.matrix(Z)

  linkfun_mu <- family$linkfun_mu
  linkinv_mu <- family$linkinv_mu
  mu.eta <- family$mu.eta
  mu2.eta2 <- family$mu2.eta2

  linkfun_phi <- family$linkfun_phi
  linkinv_phi <- family$linkinv_phi
  phi.vi <- family$phi.vi
  phi2.vi2 <- family$phi2.vi2

  variance <- family$variance

  hess_phi <- family$hess_phi
  hess_mu_phi <- family$hess_mu_phi
  loglik <- family$loglik

  lower <- c(lower_mu, lower_phi)
  upper <- c(upper_mu, upper_phi)

  map_functions <- make_map_function(lower, upper)
  map <- map_functions$map
  invert <- map_functions$invert
  jacobian <- map_functions$jacobian

  npar_mu <- length(lower_mu)
  npar_phi <- length(lower_phi)
  npar <- npar_mu + npar_phi

  sa <- sa_fit(y = y, X = X, Z = Z, family = family,
               f_mu = f_mu, J_mu = J_mu, H_mu = H_mu,
               f_phi = f_phi, J_phi = J_phi, H_phi = H_phi,
               beta_start = beta_start, lower_mu = lower_mu, upper_mu = upper_mu,
               gamma_start = gamma_start, lower_phi = lower_phi, upper_phi = upper_phi,
               mult = mult, nsim = nsim, sa_control = sa_control,
               expected = expected, verbose = verbose)

  nr <- tryCatch(gnlmsa_fit(y = y, X = X, Z = Z, family = gnlmsa_Gamma("identity"),
                            f_mu = cobb_douglas()$f, J_mu = cobb_douglas()$J,
                            H_mu = cobb_douglas()$H,
                            f_phi = Linear()$f, J_phi = Linear()$J, H_phi = Linear()$H,
                            beta_start = sa$beta, gamma_start = sa$gamma,
                            maxit = 100, tol = 1e-05, expected = TRUE),
                 error = function(e) "FAILED")

  if (!inherits(nr, "gnlmsa_fit")) {
    nr_failed <- TRUE
    fit <- sa
  } else if (nr$loglik > sa$loglik) {
    fit <- nr
    nr_better <- nr_failed <- TRUE
  } else {
    fit <- sa
    nr_better <- FALSE
    nr_failed <- FALSE
  }

  if (missing(beta_names)) {
    beta_names <- paste0("beta", 1:length(fit$beta))
  }

  if (missing(gamma_names)) {
    gamma_names <- paste0("gamma", 1:length(fit$gamma))
  }

  beta <- fit$beta
  gamma <- fit$gamma

  if (length(beta_names) != length(beta)) {
    warning("provided beta_names doesn't match the number of actual parameters, switching to default name")
    beta_names <- paste0("beta", 1:length(fit$beta))
  }

  if (length(gamma_names) != length(gamma)) {
    warning("provided gamma_names doesn't match the number of actual parameters, switching to default name")
    gamma_names <- paste0("gamma", 1:length(fit$gamma))
  }


  coef_names <- c(beta_names, gamma_names)
  names(beta) <- beta_names
  names(gamma) <- gamma_names

  loglik <- fit$loglik
  nobs <- NROW(y)
  df.residuals <- nobs - npar
  df <- npar

  out <- list(coef = c(beta, gamma),
              beta = beta,
              gamma = gamma,
              loglik = loglik,
              eta = fit$eta,
              mu = fit$mu,
              vi = fit$vi,
              phi = fit$phi,
              nsim = nsim,
              mult = mult,
              expected = expected,
              nobs = nobs,
              df = df,
              df.residuals = df.residuals,
              npar_mu = npar_mu,
              coef_names = coef_names,
              family = family,
              nr_failed = nr_failed,
              nr_better = nr_better,
              y = y,
              X = X,
              Z = Z,
              mean_model = list(f_mu = f_mu, J_mu = J_mu, H_mu = H_mu),
              dispersion_model = list(f_phi = f_phi, J_phi = J_phi, H_phi = H_phi),
              map_functions = map_functions,
              lower = lower,
              upper = upper)
  class(out) <- "gnlmsa"

  out
}
