


#-------------------------------------------------------------------------------



#' Options for Simulated Annealing algorithm
#'
#' @param iterations integer indicating the number of iterations to perform.
#' @param compute_v how many iterations wait in order to recompute variance-covariance matrix.
#' @param initial_temperature initial temperature of the algorithm.
#' @param final_temperature final temperature of the algorithm.
#' @param restart_if_stuck how many iterations wait before restarting from last best parameters value if the algorithm stucks.
#' @param save_history whether to save full history of the algorithm.
#'
#' @return a list containing the specified options.
#' @export
sa_control <- function(iterations = 1000,
                       compute_v = floor(iterations/10),
                       initial_temperature = 100,
                       final_temperature = 1,
                       restart_if_stuck = floor(iterations*0.33),
                       save_history = FALSE
) {

  if (compute_v > iterations) {
    warning("compute_v must be less or equal to iterations, setting compute_v equal to iterations")
    compute_v <- iterations
  }

  if (initial_temperature < final_temperature) {
    stop("initial temperature must be greater than final temperature")
  }

  list(iterations = iterations,
       compute_v = compute_v,
       initial_temperature = initial_temperature,
       final_temperature = final_temperature,
       restart_if_stuck = restart_if_stuck,
       save_history = save_history)
}



#-------------------------------------------------------------------------------



#' Simulated Annealing for Generalized Non-Linear Models
#'
#' @description
#' Optimizes parameters of a generalized linear model with separate functions for mean and dispersion
#' using simulated annealing. The algorithm explores the parameter space to find a candidate maximum point
#' of the log-likelihood function.
#'
#' @param y numerical vector of response values.
#' @param X design matrix for the mean model.
#' @param Z design matrix for the dispersion model.
#' @param family a `family_gnlmsa` object.
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
#' @param expected logical indicating whether to use expected (`expected = TRUE`, default) or observed (`expected = FALSE`) hessian.
#' @param verbose logical indicating whether to print progress information.
#'
#' @details
#' The function implements a Simulated Annealing algorithm for finding the values of the parameters that possibly maximize the
#' log-likelihood of a Generalized Non-Linear model. Those values can be used as starting value in order to converge to a possibly global maximum via other optimizations algorithms.
#'
#' The algorithm maps between constrained and unconstrained parameter spaces to respect the bounds
#' specified by `lower_mu`, `upper_mu`, `lower_phi`, and `upper_phi`.
#'
#' At regular intervals (controlled by `sa_control$compute_v`), the algorithm recomputes the
#' variance-covariance matrix used for the proposal distribution based on the Hessian of the
#' log-likelihood.
#'
#' @return A list with class "sa" containing:
#' \item{par}{Estimated parameters}
#' \item{beta}{Estimated parameters for the mean component}
#' \item{gamma}{Estimated parameters for the dispersion component}
#' \item{loglik}{Maximum log-likelihood value found}
#' \item{eta}{Predictor for the mean component}
#' \item{mu}{Fitted mean values}
#' \item{vi}{Predictor for the dispersion component}
#' \item{phi}{Fitted dispersion values}
#' \item{iterations}{Number of iterations performed}
#' \item{history}{Vector of log-likelihood values at each iteration (if `save_history = TRUE`)}
#' \item{control}{Control parameters used for the algorithm.}
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
#'                mult = mult, nsim = nsim, sa_control = sa_control_list,
#'                expected = expected, verbose = verbose)
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
                    mult, nsim, sa_control = sa_control(),
                    expected = TRUE, verbose = TRUE) {

  if(length(lower_mu) != length(upper_mu)) stop("lower_mu and upper_mu must have the same length")
  if(length(lower_phi) != length(upper_phi)) stop("lower_phi and upper_phi must have the same length")

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


  sa_iterations <- sa_control$iterations
  sa_compute_v <- sa_control$compute_v
  iter_compute_v <- which(seq_len(sa_iterations) %% sa_compute_v == 0)
  restart_if_stuck <- sa_control$restart_if_stuck
  save_history <- sa_control$save_history


  par0_con <- c(beta_start, gamma_start)
  par0_unc <- map(par0_con)
  beta0 <- par0_con[1:npar_mu]
  gamma0 <- par0_con[(npar_mu + 1):npar]

  eta0 <- f_mu(X, beta0)
  mu0 <- linkinv_mu(eta0)

  vi0 <- f_phi(Z, gamma0)
  phi0 <- linkinv_phi(vi0)

  l0 <- sum(loglik(y, mu0, phi0))

  v <- diag(1e-03, npar)
  j <- diag(jacobian(par0_con))
  v <- j%*%v%*%t(j)


  l_best <- l0
  par_best <- par0_con

  initial_temperature <- sa_control$initial_temperature
  final_temperature <- sa_control$final_temperature
  const_temp <- log(initial_temperature/final_temperature)/sa_iterations
  temp <- initial_temperature*exp(-const_temp*0)



  if (save_history) {
    loglik_history <- numeric(sa_iterations + 1)
    loglik_history[1] <- l0
  }

  stuck_it <- 0

  for(i in seq_len(sa_iterations)){
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

      v <- tryCatch(solve(-hess(h_mu, h_phi, h_mu_phi)), error = function(e) diag(1e-06, npar))

      if (any(diag(v) < 0) | any(is.nan(v)) | any(is.na(v))) {
        v <- diag(1e-06, npar)
      }

      j <- diag(jacobian(par0_con))
      v <- j%*%v%*%t(j)

    }


    par1_unc <- sample_par(par0_unc, v, mult, npar)
    par1_con <- invert(par1_unc)
    beta1 <- par1_con[1:npar_mu]
    gamma1 <- par1_con[(npar_mu + 1):npar]

    eta1 <- f_mu(X, beta1)
    mu1 <- linkinv_mu(eta1)

    vi1 <- f_phi(Z, gamma1)
    phi1 <- linkinv_phi(vi1)

    l1 <- sum(loglik(y, mu1, phi1))
    delta <- l1 - l0

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
  eta_best <- f_mu(X, beta_best)
  vi_best <- f_phi(Z, gamma_best)
  mu_best <- linkinv_mu(eta_best)
  phi_best <- linkinv_phi(vi_best)

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
              control = sa_control_list)

  if (save_history) {
    out$history <- loglik_history
  }

  class(out) <- "sa"

  out

}



#-------------------------------------------------------------------------------


