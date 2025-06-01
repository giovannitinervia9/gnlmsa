


#-------------------------------------------------------------------------------



#' Fit a Generalized Non-Linear Model with Mean and Dispersion Components
#'
#' Fits a Generalized Non-Linear Model (GNLM) by combining global optimization via Simulated Annealing
#' and local refinement via Newton–Raphson. The model includes both a mean component and a dispersion
#' component, each of which can be specified using user-defined nonlinear predictors.
#'
#' @description
#' The fitting procedure consists of two stages:
#' \enumerate{
#'   \item Simulated Annealing is used to explore the parameter space and identify a good candidate solution.
#'   \item Newton–Raphson is then used to refine this candidate. If the local optimization fails, or
#'         the Simulated Annealing solution yields a higher log-likelihood, the SA result is retained.
#' }
#'
#' @param y Numeric vector of response values.
#' @param mean_model A list that includes the nonlinear predictor function, Jacobian, Hessian, lower and upper bounds,
#' and design matrix for the mean component. The list should have the following components:
#' - `f`: Function for the mean component nonlinear predictor.
#' - `J`: Jacobian of the f function.
#' - `H`: Hessian of the f function.
#' - `lower`: Lower bounds for the mean component parameters.
#' - `upper`: Upper bounds for the mean component parameters.
#' - `X`: Design matrix for the mean component.
#' @param beta_start Initial values for the parameters of the mean component.
#' @param family A family_gnlmsa object specifying the distributional assumptions and required derivatives.
#' @param dispersion_model (Optional) A list that includes the nonlinear predictor function, Jacobian, Hessian,
#' lower and upper bounds, and design matrix for the dispersion component. The list should have the
#' following components:
#' - `f`: Function for the dispersion component nonlinear predictor.
#' - `J`: Jacobian of the f function.
#' - `H`: Hessian of the f function.
#' - `lower`: Lower bounds for the dispersion component parameters.
#' - `upper`: Upper bounds for the dispersion component parameters.
#' - `X`: Design matrix for the dispersion component.
#' @param gamma_start Initial values for the parameters of the dispersion component.
#' @param mult Proposal scaling factor for the Simulated Annealing algorithm.
#' @param nsim Number of independent Simulated Annealing simulations (currently only 1 is supported).
#' @param sa_control_params A list of Simulated Annealing options created via [sa_control()].
#' @param nr_control_params A list of Newton-Raphson options created via [nr_control()].
#' @param fixed_params A list containing two numeric vectors of the same length:
#'   - The first vector specifies the indices of the parameters to be fixed.
#'   - The second vector provides the corresponding fixed values for those parameters.
#' @param beta_names (Optional) Names of the parameters in the mean component.
#' @param gamma_names (Optional) Names of the parameters in the dispersion component.
#'
#' @return A list of class `"gnlmsa"` containing:
#' \describe{
#'   \item{coef}{Named vector of estimated parameters (concatenation of \code{beta} and \code{gamma}).}
#'   \item{beta, gamma}{Estimated coefficients for mean and dispersion components.}
#'   \item{loglik}{Log-likelihood of the final fitted model.}
#'   \item{eta, mu}{Predictor and fitted values for the mean component.}
#'   \item{vi, phi}{Predictor and fitted values for the dispersion component.}
#'   \item{nsim, mult}{Settings used for Simulated Annealing.}
#'   \item{nobs, df, df.residuals}{Sample size and degrees of freedom.}
#'   \item{npar_mu}{Number of parameters for the mean component.}
#'   \item{nr_failed}{Logical; \code{TRUE} if the Newton–Raphson step failed.}
#'   \item{nr_better}{Logical; \code{TRUE} if Newton–Raphson log-likelihood was better than Simulated Annealing.}
#'   \item{mean_model, dispersion_model}{List of functions used to specify model components.}
#'   \item{map_functions}{List of mapping functions used to transform bounded to unbounded parameters.}
#'   \item{lower, upper}{Parameter bounds.}
#'   \item{family}{The `family_gnlmsa` object used.}
#'   \item{X, Z, y}{Original data matrices and response vector.}
#'   \item{one_parameter}{Logical indicating wether the family of the object is a one parameter family.}
#'   \item{gradient, hessian}{Gradient and Hessian at the estimated values of the parameters (free
#'         parameters only).}
#'   \item{nr}{The object returned from `gnlmsa_fit()` function.}
#'   \item{sa}{The object returned from `sa_fit()` function.}
#' }
#' @export
#' @importFrom stats var
#' @examples
#' \dontrun{
#' y <- productivity$Y
#' X <- productivity[, -1]
#' beta_start <- c(5, .5, .5)
#' fit <- gnlmsa(y, mean_model = cobb_douglas(X = X), beta_start, family = gnlmsa_Gamma("identity"),
#'               sa_control_params = sa_control(1000))
#' }
gnlmsa <- function (y, mean_model, beta_start,
                    family = gnlmsa_gaussian(), dispersion_model = NULL, gamma_start,
                    mult = 1, nsim = 1, sa_control_params = sa_control(),
                    nr_control_params = nr_control(),
                    fixed_params = NULL,
                    beta_names, gamma_names) {

  one_parameter <- family$family %in% c("gnlmsa_poisson", "gnlmsa_binomial")

  nobs <- NROW(y)

  f_mu <- mean_model$f
  J_mu <- mean_model$J
  H_mu <- mean_model$H
  lower_mu <- mean_model$lower
  upper_mu <- mean_model$upper
  X <- as.matrix(mean_model$X)

  if (is.null(dispersion_model) | one_parameter) {
    f_phi <- Linear()$f
    J_phi <- Linear()$J
    H_phi <- Linear()$H
    lower_phi <- -Inf
    upper_phi <- Inf
    Z <- as.matrix(rep(1, nobs))
    if (missing(gamma_start)) gamma_start <- log(var(y))
  } else {
    f_phi <- dispersion_model$f
    J_phi <- dispersion_model$J
    H_phi <- dispersion_model$H
    lower_phi <- dispersion_model$lower
    upper_phi <- dispersion_model$upper
    Z <- dispersion_model$X
  }

  if(length(lower_mu) != length(upper_mu)) stop("lower_mu and upper_mu must have the same length")
  if(length(lower_phi) != length(upper_phi)) stop("lower_phi and upper_phi must have the same length")

  if(any(beta_start <= lower_mu) | any(beta_start >= upper_mu)) {
    stop("Starting values for beta are outside the lower and upper bounds defined.")
  }

  if(any(gamma_start <= lower_phi) | any(gamma_start >= upper_phi)) {
    stop("Starting values for gamma are outside the lower and upper bounds defined.")
  }

  if(is.null(J_mu)) J_mu <- make_jacobian(f_mu)
  if(is.null(H_mu)) H_mu <- make_hessian(f_mu)
  if(is.null(J_phi)) J_phi <- make_jacobian(f_phi)
  if(is.null(H_phi)) H_phi <- make_hessian(f_phi)

  lower <- c(lower_mu, lower_phi)
  upper <- c(upper_mu, upper_phi)

  npar_mu <- length(lower_mu)
  npar_phi <- length(lower_phi)
  npar <- npar_mu + npar_phi


  sa <- sa_fit(y = y, X = X, Z = Z, family = family,
               f_mu = f_mu, J_mu = J_mu, H_mu = H_mu,
               f_phi = f_phi, J_phi = J_phi, H_phi = H_phi,
               beta_start = beta_start, lower_mu = lower_mu, upper_mu = upper_mu,
               gamma_start = gamma_start, lower_phi = lower_phi, upper_phi = upper_phi,
               fixed_params = fixed_params,
               mult = mult, nsim = nsim, sa_control_params = sa_control_params)

  nr <- tryCatch(gnlmsa_fit(y = y, X = X, Z = Z, family = family,
                            f_mu = f_mu, J_mu = J_mu,
                            H_mu = H_mu,
                            f_phi = f_phi, J_phi = J_phi, H_phi = H_phi,
                            beta_start = sa$beta, lower_mu = lower_mu, upper_mu = upper_mu,
                            gamma_start = sa$gamma, lower_phi = lower_phi, upper_phi = upper_phi,
                            fixed_params = fixed_params,
                            nr_control_params = nr_control_params),
                 error = function(e) "FAILED")

  if (!inherits(nr, "gnlmsa_fit")) {
    nr_failed <- TRUE
    nr_better <- FALSE
    fit <- sa
    warning("\nNewton-Raphson optimization failed to optimize log-likelihood function.\n Try to use more iterations for Simulated Annealing algorithm.")
  } else if (nr$loglik >= sa$loglik) {
    fit <- nr
    nr_better <- TRUE
    nr_failed <- FALSE
  } else {
    fit <- sa
    nr_better <- FALSE
    nr_failed <- FALSE
    warning("\nNewton-Raphson optimization failed to optimize log-likelihood function.\n Try to use more iterations for Simulated Annealing algorithm.")
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

  if (one_parameter) {
    df <- length(beta)
  } else {
    df <- npar
  }

  df.residuals <- nobs - df


  out <- list(coefficients = c(beta, gamma),
              beta = beta,
              gamma = gamma,
              loglik = loglik,
              eta = fit$eta,
              mu = fit$mu,
              vi = fit$vi,
              phi = fit$phi,
              nsim = nsim,
              mult = mult,
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
              map_functions = sa$map_functions,
              lower = lower,
              upper = upper,
              fixed_params = fixed_params,
              one_parameter = one_parameter,
              gradient = fit$gradient,
              hessian = fit$hessian,
              nr = nr,
              sa = sa)
  class(out) <- "gnlmsa"

  out
}
