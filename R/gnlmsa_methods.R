


#-------------------------------------------------------------------------------



#' Print Method for `gnlmsa` objects
#'
#' Provides a customized printout for objects of class `"gnlmsa"`, summarizing the fitted
#' model information, the estimated coefficients for both the mean and dispersion components,
#' and the status of the Newton–Raphson optimization.
#'
#' @param x An object of class `"gnlmsa"`, typically returned by [gnlmsa()].
#' @param digits Integer. The number of significant digits to use when printing.
#'   Defaults to \code{max(3L, getOption("digits") - 3L)}.
#' @param ... Further arguments passed to or from other methods (currently unused).
#'
#' @details
#' The output includes:
#' \itemize{
#'   \item The model type and family name.
#'   \item The estimated coefficients for the mean component.
#'   \item The estimated coefficients for the dispersion component.
#'   \item An indication of whether the Newton–Raphson optimization step succeeded or failed.
#' }
#' If Newton–Raphson optimization fails, the function reports the estimates obtained via Simulated Annealing instead.
#'
#' @return The object \code{x} is returned invisibly.
#'
#'
#' @export
print.gnlmsa <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  cat("\nGeneralized Non-Linear Model")
  cat(paste0("\nFamily: ", x$family$family), "\n")
  cat("\nCoefficients for mean component:\n")
  print.default(format(x$beta), digits = digits, print.gap = 2L, quote = FALSE)
  cat("\nCoefficients for dispersion component:\n")
  print.default(format(x$gamma), digits = digits, print.gap = 2L, quote = FALSE)

  if (!x$nr_failed) {
    cat("\nNewton-Raphson optimization: SUCCEEDED\nReporting Newton-Raphson estimates")
  } else {
    cat("\nNewton-Raphson optimization: FAILED\nReporting Simulated Annealing estimates")
  }
  invisible(x)

}



#-------------------------------------------------------------------------------



#' Variance-Covariance Matrix for Generalized Non-Linear Model Objects
#'
#' Computes the estimated variance-covariance matrix of the parameter estimates for a fitted
#' Generalized Non-Linear Model (GNLM) object of class `"gnlmsa"`.
#'
#' @param object An object of class `"gnlmsa"`, typically returned by [gnlmsa()].
#' @param expected Logical. If \code{TRUE} (default), the expected Fisher Information matrix is used
#'   to compute the variance-covariance matrix. If \code{FALSE}, the observed Fisher Information matrix
#'   is used instead.
#' @param ... Further arguments passed to or from other methods (currently unused).
#'
#' @details
#' The function computes the Hessian matrix of the log-likelihood with respect to both mean and dispersion parameters,
#' optionally using either the expected or the observed Fisher Information or the observed Hessian.
#'
#' If the Fisher Information matrix is not invertible, a matrix filled with \code{NA} values is returned instead, and no error is thrown.
#'
#' Column and row names of the resulting matrix correspond to the names of the model parameters.
#'
#' @return A symmetric numeric matrix representing the estimated variance-covariance matrix
#'   of the parameter estimates.
#'
#' @export
vcov.gnlmsa <- function(object, expected = TRUE, ...) {

  y <- object$y

  X <- object$X
  beta <- object$beta
  mu <- object$mu
  eta <- object$eta
  f_mu <- object$mean_model$f_mu
  J_mu <- object$mean_model$J_mu
  H_mu <- object$mean_model$H_mu


  Z <- object$Z
  gamma <- object$gamma
  phi <- object$phi
  vi <- object$vi
  f_phi <- object$dispersion_model$f_phi
  J_phi <- object$dispersion_model$J_phi
  H_phi <- object$dispersion_model$H_phi

  family <- object$family
  mu.eta <- family$mu.eta
  mu2.eta2 <- family$mu2.eta2
  phi.vi <- family$phi.vi
  phi2.vi2 <- family$phi2.vi2
  variance <- family$variance

  hess_phi <- family$hess_phi
  hess_mu_phi <- family$hess_mu_phi


  h_mu <- hess_mu(y, X, beta, mu, eta, phi, f_mu, J_mu, H_mu,
                  mu.eta, mu2.eta2, variance, expected)

  h_phi <- hess_phi(y, Z, gamma, phi, vi, mu, f_phi, J_phi, H_phi,
                    phi.vi, phi2.vi2, expected)

  h_mu_phi <- hess_mu_phi(y, X, Z, beta, gamma, mu, eta, phi, vi, f_mu,
                          J_mu, f_phi, J_phi, mu.eta, phi.vi, expected)

  v <- tryCatch(solve(-hess(h_mu, h_phi, h_mu_phi)),
                error = function(e) {
                  matrix(NA, length(c(beta, gamma)), length(c(beta, gamma)))
                }
  )

  colnames(v) <- rownames(v) <- object$coef_names

  v

}



#-------------------------------------------------------------------------------



#' Confidence Intervals for Generalized Non-Linear Model Parameters
#'
#' Computes confidence intervals for the parameters of a fitted Generalized Non-Linear Model (GNLM) object
#' of class `"gnlmsa"`, based on the variance-covariance matrix of the parameter estimates.
#'
#' @param object An object of class `"gnlmsa"` returned by [gnlmsa()].
#' @param parm Parameters for which confidence intervals are required. Currently ignored; confidence intervals are computed for all parameters.
#' @param level Confidence level required. Default is \code{0.95} for 95% confidence intervals.
#' @param test Character string specifying the type of test to use. Currently, only `"Wald"` is supported.
#'   If a different test is specified, a warning is issued and the function falls back to `"Wald"`.
#' @param expected Logical. If \code{TRUE} (default), uses the expected Fisher information to compute standard errors;
#'   otherwise uses the observed information.
#' @param ... Further arguments passed to or from other methods (currently unused).
#'
#' @details
#' Currently, only Wald-type confidence intervals are implemented.
#'
#' Confidence intervals are computed on the transformed (unconstrained) scale and then mapped back
#' to the original (possibly constrained) parameter space.
#'
#' The variance-covariance matrix used in the computation accounts for the parameter transformations,
#' via the Jacobian of the transformation at the point estimates.
#'
#' @return A data frame with two columns, giving the lower and upper bounds of the confidence intervals
#' for each parameter.
#'
#'
#' @export
confint.gnlmsa <- function(object, parm, level = 0.95, test = c("Wald", "Rao", "LRT"), expected = TRUE, ...) {

  test <- match.arg(test)

  if (test != "Wald") {
    test <- "Wald"
    warning("Currently only Wald confindence intervals are supported. Switching to Wald confidence intervals.")
  }

  if (test == "Wald") {

    par <- object$coefficients
    v <- vcov(object, expected)
    map_functions <- object$map_functions

    map <- map_functions$map
    invert <- map_functions$invert
    map_jacobian <- map_functions$map_jacobian
    map_functions

    j <- map_jacobian(par)

    v_map <- diag(j)%*%v%*%t(diag(j))
    se_map <- sqrt(diag(v_map))
    sig.level <- 1 - level
    q <- qnorm(1 - sig.level/2)

    inf <- invert(map(par) - q*se_map)
    sup <- invert(map(par) + q*se_map)


    out <- data.frame(inf = inf, sup = sup)
    colnames(out) <- paste0(colnames(out), "_", paste0(c(sig.level/2*100, (1 - sig.level/2)*100), "%"))
    rownames(out) <- object$coef_names


  }

  out

}



#-------------------------------------------------------------------------------



#' Extract Model Coefficients from a Generalized Non-Linear Model
#'
#' Returns the estimated model coefficients from a fitted Generalized Non-Linear Model (GNLM)
#' object of class `"gnlmsa"`.
#'
#' @param object An object of class `"gnlmsa"`, typically returned by [gnlmsa()].
#' @param ... Further arguments passed to or from other methods (currently unused).
#'
#' @details
#' This function extracts and returns the estimated coefficients for both the mean and dispersion components,
#' combined into a single named numeric vector.
#'
#' The coefficients are returned on the original (constrained) parameter space.
#'
#' @return A named numeric vector of model coefficients.
#'
#' @export
coef.gnlmsa <- function(object, ...) {
  object$coefficients
}



#-------------------------------------------------------------------------------

















