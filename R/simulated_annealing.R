


#-------------------------------------------------------------------------------



#' Control Parameters for the Simulated Annealing Algorithm
#'
#' Defines a set of user-configurable options to control the behavior of the Simulated Annealing (SA)
#' routine used in `gnlmsa` for global exploration of the parameter space in Generalized Non-Linear Models.
#'
#' These settings govern the temperature schedule, frequency of variance updates, restart conditions,
#' and tracking behavior during the annealing process.
#'
#' @param iterations Integer. Total number of iterations for the Simulated Annealing algorithm.
#' @param initial_temperature Numeric. Initial value of the temperature. Should be strictly greater than \code{final_temperature}.
#' @param final_temperature Numeric. Final value of the temperature at which the cooling schedule stops.
#' @param chol Logical. If \code{TRUE}, sampling is performed using the Cholesky decomposition of v, otherwhise `mvtnorm::rmvnorm()` is used.
#' @param update_v Logical. If \code{TRUE}, the variance covariance matrix of the sampler is updated. Default is \code{FALSE}.
#' @param compute_v Integer. Number of iterations after which the variance-covariance matrix used
#'   for proposal generation is updated. Default is every 30% of total iterations.
#' @param restart_if_stuck Integer. Number of consecutive iterations without improvement after which
#'   the algorithm restarts from the last best parameter vector. Helps escape local optima.
#' @param save_history Logical. If \code{TRUE}, the full parameter and log-likelihood trajectory is stored.
#' @param expected Logical; if \code{TRUE}, the expected (Fisher) information matrix is used to update the proposal variance.
#' @param verbose Logical; if \code{TRUE}, prints algorithm progress every 100 iterations.
#' @return A list of class `"sa_control"` containing the control parameters to be used in the SA routine.
#'
#' @details
#' These parameters influence the convergence and robustness of the SA-based global search.
#' After reaching the final temperature or maximum iteration count, the best solution found
#' can be passed to a Newton-Raphson or other local optimization method for refinement.
#'
#' @export
sa_control <- function(iterations = 1000,
                       initial_temperature = 100,
                       final_temperature = 1,
                       chol = T,
                       update_v = FALSE,
                       compute_v = floor(iterations*0.3),
                       restart_if_stuck = floor(iterations*0.3),
                       save_history = FALSE,
                       expected = TRUE,
                       verbose = TRUE
) {

  if (compute_v > iterations) {
    warning("compute_v must be less or equal to iterations, setting compute_v equal to iterations")
    compute_v <- iterations
  }

  if (initial_temperature < final_temperature) {
    stop("initial temperature must be greater than final temperature")
  }

  list(iterations = iterations,
       initial_temperature = initial_temperature,
       final_temperature = final_temperature,
       chol = chol,
       update_v = update_v,
       compute_v = compute_v,
       restart_if_stuck = restart_if_stuck,
       save_history = save_history,
       expected = expected,
       verbose = verbose)
}



#-------------------------------------------------------------------------------



#' Simulated Annealing Optimization for Generalized Non-Linear Models
#'
#' Performs parameter estimation for Generalized Non-Linear Models (GNLMs) using a Simulated Annealing (SA) algorithm.
#' This global optimization technique searches for a maximum of the log-likelihood function by exploring the
#' parameter space in a probabilistic way, before passing control to a local optimizer.
#'
#' @param y Numeric vector of response values.
#' @param X Numeric design matrix for the mean component.
#' @param Z Numeric design matrix for the dispersion component.
#' @param family An object of class [`family_gnlmsa`] that defines the distribution, link functions, and relevant derivatives.
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
#' @param fixed_params A list containing two numeric vectors of the same length:
#'   - The first vector specifies the indices of the parameters to be fixed.
#'   - The second vector provides the corresponding fixed values for those parameters.
#' @param mult Multipliers for the proposal variance in the SA algorithm.
#' @param nsim Number of independent SA chains to run (for potential parallelization; currently only 1 is supported).
#' @param sa_control_params A list of control parameters for the SA routine, as returned by [sa_control()].
#'
#' @details
#' The SA algorithm operates on a transformed parameter space: all constrained parameters are mapped
#' to \eqn{\mathbb{R}} via bijective transformations to ensure feasibility.
#'
#' At intervals determined by `sa_control_params$compute_v`, the algorithm estimates a local variance-covariance
#' matrix from the (expected or observed) Hessian of the log-likelihood and uses it to adapt the proposal distribution.
#'
#' If no improvement is detected after a certain number of iterations (`restart_if_stuck`), the algorithm
#' restarts from the best parameter configuration found so far. At the end of the annealing schedule,
#' the solution can be passed to a local optimizer (e.g., Newton-Raphson).
#'
#' @return A list of class `"sa"` with elements:
#' \describe{
#'   \item{par}{Estimated parameter vector (concatenated \code{beta} and \code{gamma}).}
#'   \item{beta, gamma}{Estimated parameters for the mean and dispersion components, respectively.}
#'   \item{loglik}{Final log-likelihood value.}
#'   \item{eta, mu}{Nonlinear predictor and fitted values for the mean component.}
#'   \item{vi, phi}{Nonlinear predictor and fitted values for the dispersion component.}
#'   \item{iterations}{Total number of iterations performed.}
#'   \item{history}{Optional. Vector of log-likelihood values at each iteration, if \code{save_history = TRUE}.}
#'   \item{control}{Control parameters used during optimization.}
#'   \item{map_functions}{List of mapping functions used to transform bounded to unbounded parameters.}
#'   \item{gradient, hessian}{Gradient and Hessian at the estimated values of the parameters (free
#'         parameters only).}
#' }
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
#' beta <- c(3, .4, .4)
#' eta <- f_mu(X, beta)
#' mu <- fam_gamma$linkinv_mu(eta)
#'
#' # parameters for dispersion component
#' gamma <- c(0.5, .01, .05)
#' vi <- f_phi(Z, gamma)
#' phi <- fam_gamma$linkinv_phi(vi)
#'
#' # simulate response
#' y <- rgamma(n, phi, phi/mu)
#' try(pairs(cbind(X, y)))
#' summary(cbind(X, y))
#'
#' # arguments
#' y = y # y
#' X = X # X
#' Z = Z # Z
#' family = gnlmsa_Gamma(link_mu = "identity", link_phi = "log") # family
#' f_mu = cobb_douglas()$f
#' J_mu = cobb_douglas()$J
#' H_mu = cobb_douglas()$H
#' f_phi = Linear()$f
#' J_phi = Linear()$J
#' H_phi = Linear()$H
#' beta_start = c(1, .9, .8)
#' lower_mu = c(0, 0, 0) # lower_mu
#' upper_mu = c(Inf, Inf, Inf) # upper_mu
#' gamma_start = c(1, 0, 0)
#' lower_phi = c(-Inf, -Inf, -Inf) # lower_phi
#' upper_phi = c(Inf, Inf, Inf) # upper_phi
#' mult = 1 # mult
#' nsim = 1 # nsim (number of SA simulations to run in parallel)
#' sa_control_list = sa_control(iterations = 1000,
#'                              initial_temperature = 100,
#'                              final_temperature = .01,
#'                              restart_if_stuck = floor(1000 * 0.01)) # sa_control
#' expected = TRUE
#' verbose = TRUE
#'
#'
#' fit1 <- sa_fit(y = y, X = X, Z = Z, family = family,
#'                f_mu = f_mu, J_mu = J_mu, H_mu = H_mu,
#'                f_phi = f_phi, J_phi = J_phi, H_phi = H_phi,
#'                beta_start = beta_start, lower_mu = lower_mu, upper_mu = upper_mu,
#'                gamma_start = gamma_start, lower_phi = lower_phi, upper_phi = upper_phi,
#'                mult = mult, nsim = nsim, sa_control_params = sa_control_list)
#'
#' fit1
#' }
#'
#' @seealso [sa_control()] for creating control parameters
#'
#'
#' @importFrom stats runif
#' @export
sa_fit <- function (y, X, Z, family,
                    f_mu, J_mu, H_mu,
                    f_phi, J_phi, H_phi,
                    beta_start, lower_mu, upper_mu,
                    gamma_start, lower_phi, upper_phi,
                    fixed_params = NULL,
                    mult, nsim, sa_control_params = sa_control()) {

  if(length(lower_mu) != length(upper_mu)) stop("lower_mu and upper_mu must have the same length")
  if(length(lower_phi) != length(upper_phi)) stop("lower_phi and upper_phi must have the same length")

  if(missing(J_mu)) J_mu <- make_jacobian(f_mu)
  if(missing(H_mu)) H_mu <- make_hessian(f_mu)
  if(missing(J_phi)) J_phi <- make_jacobian(f_phi)
  if(missing(H_phi)) H_phi <- make_hessian(f_phi)

  one_parameter <- family$family %in% c("gnlmsa_poisson", "gnlmsa_binomial")

  npar_mu <- length(lower_mu)
  npar_phi <- length(lower_phi)
  npar <- npar_mu + npar_phi

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

  loglik <- family$loglik
  grad_mu <- family$grad_mu
  hess_mu <- family$hess_mu
  grad_phi <- family$grad_phi
  hess_phi <- family$hess_phi
  hess_mu_phi <- family$hess_mu_phi

  lower <- c(lower_mu, lower_phi)[free_par]
  upper <- c(upper_mu, upper_phi)[free_par]

  map_functions <- make_map_function(lower, upper)
  map <- map_functions$map
  invert <- map_functions$invert
  jacobian <- map_functions$map_jacobian




  sa_iterations <- sa_control_params$iterations
  use_chol <- sa_control_params$chol
  sa_compute_v <- sa_control_params$compute_v
  sa_update_v <- sa_control_params$update_v
  iter_compute_v <- which(seq_len(sa_iterations) %% sa_compute_v == 0)
  restart_if_stuck <- sa_control_params$restart_if_stuck
  save_history <- sa_control_params$save_history
  verbose <- sa_control_params$verbose
  expected <- sa_control_params$expected

  par0_con <- c(beta_start, gamma_start)
  par0_unc <- par1_con <- par1_unc <- numeric(npar)

  if (fixed) {
    par0_con[fixed_position] <- fixed_values
    par0_unc[fixed_position] <- fixed_values
    par1_con[fixed_position] <- fixed_values
    par1_unc[fixed_position] <- fixed_values
  }


  par0_unc[free_par] <- map(par0_con[free_par])
  beta0 <- par0_con[1:npar_mu]
  gamma0 <- par0_con[(npar_mu + 1):npar]

  eta0 <- drop(f_mu(X, beta0))
  mu0 <- linkinv_mu(eta0)

  vi0 <- drop(f_phi(Z, gamma0))
  phi0 <- linkinv_phi(vi0)

  l0 <- sum(loglik(y, mu0, phi0))

  v <- diag(1e-03, npar_free)
  j <- diag(jacobian(par0_con))[free_par, free_par]
  v <- j%*%v%*%t(j)
  if (use_chol) {
    L <- chol(v)
  }



  l_best <- l0
  par_best <- par0_con

  initial_temperature <- sa_control_params$initial_temperature
  final_temperature <- sa_control_params$final_temperature
  const_temp <- log(initial_temperature/final_temperature)/sa_iterations
  temp <- initial_temperature*exp(-const_temp*0)



  if (save_history) {
    loglik_history <- numeric(sa_iterations + 1)
    loglik_history[1] <- l0
  }

  stuck_it <- 0

  if (verbose) {
    cat("\n")
    cat("Simulated annealing:\n")
  }


  for(i in seq_len(sa_iterations)){

    if (sa_update_v) {
      # recompute variance covariance matrix
      if (i %in% iter_compute_v) {
        h_mu <- hess_mu(y, X, beta0, mu0, eta0, phi0,
                        f_mu, J_mu, H_mu,
                        mu.eta, mu2.eta2, variance, expected)

        h_phi <- hess_phi(y, Z, gamma0, phi0, vi0, mu0,
                          f_phi, J_phi, H_phi,
                          phi.vi, phi2.vi2, expected)

        h_mu_phi <- hess_mu_phi(y, X, Z, beta0, gamma0,
                                mu0, eta0, phi0, vi0,
                                f_mu, J_mu, f_phi, J_phi,
                                mu.eta, phi.vi, expected)

        h <- hess(h_mu, h_phi, h_mu_phi)[free_par, free_par]

        v <- tryCatch(solve(-h), error = function(e) diag(1e-06, npar))

        if (any(diag(v) < 0) | any(is.nan(v)) | any(is.na(v))) {
          v <- diag(1e-06, npar_free)
        }

        j <- diag(jacobian(par0_con))[free_par, free_par]
        v <- j%*%v%*%t(j)
        if (use_chol) {
          L <- chol(v)
        }
      }
    }

    if (use_chol) {
      par1_unc[free_par] <- par0_unc[free_par] + mult*drop(L%*%rnorm(npar_free))
    } else {
      par1_unc[free_par] <- sample_par(par0_unc[free_par], v, mult, npar_free)
    }
    par1_con[free_par] <- invert(par1_unc[free_par])

    if (fixed) {
      par1_con[fixed_position] <- fixed_values
    }

    beta1 <- par1_con[1:npar_mu]
    gamma1 <- par1_con[(npar_mu + 1):npar]

    eta1 <- drop(f_mu(X, beta1))
    mu1 <- linkinv_mu(eta1)

    vi1 <- drop(f_phi(Z, gamma1))
    phi1 <- linkinv_phi(vi1)

    l1 <- sum(loglik(y, mu1, phi1))
    delta <- l1 - l0

    if (is.nan(delta)) {
      iter_compute_v <- sa_iterations + 100
      v <- diag(1e-06, npar_free)
      next
    }

    if (delta > 0) {
      par0_unc <- par1_unc
      l0 <- l1
      stuck_it <- 0
    } else {
      temp <- initial_temperature*exp(-const_temp*i)
      p <- exp(delta/temp)
      if (runif(1) < p) {
        par0_unc <- par1_unc
        l0 <- l1
        stuck_it <- 0
      } else {
        stuck_it <- stuck_it + 1
      }
    }

    if (l0 > l_best) {
      par_best <- par1_con
      l_best <- l0
    }

    if (save_history) {
      loglik_history[i + 1] <- l0
    }

    if (stuck_it == restart_if_stuck) {
      par0_unc <- map(par_best)
      l0 <- l_best
    }


    if (verbose && i %% 100 == 0) {
      cat(sprintf("Iter %d, LogLik = %.4f, Temp = %.4f\n", i, l0, temp))
    }


  }

  beta_best <- par_best[1:npar_mu]
  gamma_best <- par_best[(npar_mu + 1):npar]
  eta_best <- drop(f_mu(X, beta_best))
  vi_best <- drop(f_phi(Z, gamma_best))
  mu_best <- drop(linkinv_mu(eta_best))
  phi_best <- drop(linkinv_phi(vi_best))

  g_mu <- grad_mu(y, X, beta_best, mu_best, eta_best, phi_best, f_mu, J_mu, mu.eta, variance)
  g_phi <- grad_phi(y, Z, gamma_best, phi_best, vi_best, mu_best, f_phi, J_phi, phi.vi)
  g <- c(g_mu, g_phi)[free_par]

  h_mu <- hess_mu(y, X, beta_best, mu_best, eta_best, phi_best, f_mu, J_mu, H_mu, mu.eta, mu2.eta2,
                  variance, expected)
  h_phi <- hess_phi(y, Z, gamma_best, phi_best, vi_best, mu_best, f_phi, J_phi, H_phi, phi.vi,
                    phi2.vi2, expected)
  h_mu_phi <- hess_mu_phi(y, X, Z, beta_best, gamma_best, mu_best, eta_best, phi_best, vi_best, f_mu, J_mu,
                          f_phi, J_phi, mu.eta, phi.vi, expected)
  h <- hess(h_mu, h_phi, h_mu_phi)[free_par, free_par]



  out <- list(par = par_best,
              beta = beta_best,
              gamma = gamma_best,
              loglik = l_best,
              eta = eta_best,
              mu = mu_best,
              vi = vi_best,
              phi = phi_best,
              iterations = sa_iterations,
              history = NULL,
              control = sa_control_params,
              map_functions = map_functions,
              gradient = g,
              hessian = h)

  if (save_history) {
    out$history <- loglik_history
  }

  class(out) <- "sa"

  out
}



#-------------------------------------------------------------------------------


