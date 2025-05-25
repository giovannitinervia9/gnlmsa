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

  if (!x$one_parameter) {
    cat("\nCoefficients for dispersion component:\n")
    print.default(format(x$gamma), digits = digits, print.gap = 2L, quote = FALSE)
  }

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

  one_parameter <- object$one_parameter


  hess_mu <- family$hess_mu
  hess_phi <- family$hess_phi
  hess_mu_phi <- family$hess_mu_phi


  h_mu <- hess_mu(
    y, X, beta, mu, eta, phi, f_mu, J_mu, H_mu,
    mu.eta, mu2.eta2, variance, expected
  )

  h_phi <- hess_phi(
    y, Z, gamma, phi, vi, mu, f_phi, J_phi, H_phi,
    phi.vi, phi2.vi2, expected
  )

  h_mu_phi <- hess_mu_phi(
    y, X, Z, beta, gamma, mu, eta, phi, vi, f_mu,
    J_mu, f_phi, J_phi, mu.eta, phi.vi, expected
  )

  if (one_parameter) {
    h <- h_mu
  } else {
    h <- hess(h_mu, h_phi, h_mu_phi)
  }

  v <- tryCatch(solve(-h),
    error = function(e) {
      matrix(NA, length(c(beta, gamma)), length(c(beta, gamma)))
    }
  )

  if (one_parameter) {
    colnames(v) <- rownames(v) <- names(object$beta)
  } else {
    colnames(v) <- rownames(v) <- object$coef_names
  }



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
  family <- object$family
  one_parameter <- object$one_parameter

  if (test != "Wald") {
    test <- "Wald"
    warning("Currently only Wald confindence intervals are supported. Switching to Wald confidence intervals.")
  }

  if (test == "Wald") {
    if (one_parameter) {
      par <- object$beta
      coef_names <- names(par)
    } else {
      par <- object$coefficients
      coef_names <- object$coef_names
    }

    v <- vcov.gnlmsa(object, expected)
    map_functions <- object$map_functions

    map <- map_functions$map
    invert <- map_functions$invert
    map_jacobian <- map_functions$map_jacobian

    j <- map_jacobian(par)

    v_map <- diag(j) %*% v %*% t(diag(j))
    se_map <- sqrt(diag(v_map))
    sig.level <- 1 - level
    q <- qnorm(1 - sig.level / 2)

    inf <- invert(map(par) - q * se_map)
    sup <- invert(map(par) + q * se_map)


    out <- data.frame(inf = inf, sup = sup)
    colnames(out) <- paste0(colnames(out), "_", paste0(c(sig.level / 2 * 100, (1 - sig.level / 2) * 100), "%"))
    rownames(out) <- coef_names
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
  if (object$one_parameter) {
    object$beta
  } else {
    object$coefficients
  }
}



#-------------------------------------------------------------------------------



#' Extract Fitted Values from a Generalized Non-Linear Model
#'
#' Returns fitted values from a fitted Generalized Non-Linear Model (GNLM) object of class `"gnlmsa"`,
#' including the mean component, the dispersion component, and the variance function evaluated at the fitted values.
#'
#' @param object An object of class `"gnlmsa"`, typically returned by [gnlmsa()].
#' @param type Character string specifying the type of fitted values to return.
#'   Must be one of \code{"mu"} (default), \code{"phi"}, \code{"variance"}, or \code{"all"}.
#' @param ... Further arguments passed to or from other methods (currently unused).
#'
#' @details
#' Available options for \code{type} are:
#' \itemize{
#'   \item \code{"mu"}: Fitted values for the mean component \eqn{\hat{\mu}}.
#'   \item \code{"phi"}: Fitted values for the dispersion component \eqn{\hat{\phi}}.
#'   \item \code{"variance"}: Variance evaluated at the fitted mean and dispersion, using the `variance()` function
#'         defined by the family.
#'   \item \code{"all"}: A data frame containing \code{mu}, \code{phi}, and \code{variance}.
#' }
#'
#' @return Depending on \code{type}:
#' \itemize{
#'   \item A numeric vector of fitted values (\code{"mu"}, \code{"phi"}, or \code{"variance"}).
#'   \item A data frame containing all three components if \code{type = "all"}.
#' }
#'
#'
#' @export
fitted.gnlmsa <- function(object, type = c("mu", "phi", "variance", "all"), ...) {
  type <- match.arg(type)

  if (type == "mu") {
    object$mu
  } else if (type == "phi") {
    object$phi
  } else if (type == "variance") {
    object$family$variance(object$mu, object$phi)
  } else {
    data.frame(
      mu = object$mu,
      phi = object$phi,
      variance = object$family$variance(object$mu, object$phi)
    )
  }
}



#-------------------------------------------------------------------------------



#' Extract Log-Likelihood from a Generalized Non-Linear Model
#'
#' Returns the log-likelihood of a fitted Generalized Non-Linear Model (GNLM) object of class `"gnlmsa"`.
#'
#' @param object An object of class `"gnlmsa"`, typically returned by [gnlmsa()].
#' @param ... Further arguments passed to or from other methods (currently unused).
#'
#' @details
#' The function returns the final log-likelihood value achieved during model fitting.
#' Attributes `"df"` (degrees of freedom, equal to the number of estimated parameters)
#' and `"nobs"` (number of observations) are attached to the returned object.
#'
#' The result is assigned class `"logLik"` to ensure compatibility with generic model evaluation
#' functions such as [AIC()], [BIC()], and [anova()].
#'
#' @return An object of class `"logLik"`, which is a numeric scalar representing the
#'   maximized log-likelihood, with attributes:
#'   \itemize{
#'     \item \code{df}: Degrees of freedom (number of parameters estimated).
#'     \item \code{nobs}: Number of observations used to fit the model.
#'   }
#'
#' @seealso [gnlmsa()], [AIC()], [BIC()]
#'
#' @export
logLik.gnlmsa <- function(object, ...) {
  val <- object$loglik
  df <- object$df
  attr(val, "df") <- df
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}



#-------------------------------------------------------------------------------



#' Summarize a Fitted Generalized Non-Linear Model
#'
#' Produces a summary for a fitted Generalized Non-Linear Model (GNLM) object of class `"gnlmsa"`,
#' including estimates, standard errors, confidence intervals, and model fit statistics.
#'
#' @param object An object of class `"gnlmsa"`, typically returned by [gnlmsa()].
#' @param level Confidence level for the confidence intervals (default is \code{0.95} for 95% intervals).
#' @param test Character string specifying the type of test for confidence intervals.
#'   Currently, only `"Wald"` is supported. If another test is specified, it is ignored with a warning.
#' @param expected Logical. If \code{TRUE} (default), uses the expected Fisher information to compute standard errors and intervals;
#'   otherwise uses the observed information.
#' @param ... Further arguments passed to or from other methods (currently unused).
#'
#' @details
#' The summary includes:
#' \itemize{
#'   \item Parameter estimates for the mean (\eqn{\mu}) and dispersion (\eqn{\phi}) components.
#'   \item Standard errors computed from the variance-covariance matrix.
#'   \item Confidence intervals for the parameters based on the Wald method.
#'   \item Log-likelihood, Akaike Information Criterion (AIC), and Bayesian Information Criterion (BIC).
#'   \item An indication of whether Newton–Raphson optimization succeeded or fallback to Simulated Annealing was necessary.
#' }
#'
#' @return A list of class `"summary.gnlmsa"` containing:
#' \describe{
#'   \item{coefficients.mu}{A data frame with estimates, standard errors, and confidence intervals for the mean component.}
#'   \item{coefficients.phi}{A data frame with estimates, standard errors, and confidence intervals for the dispersion component.}
#'   \item{nr_failed}{Logical; \code{TRUE} if Newton–Raphson optimization failed.}
#'   \item{nr_better}{Logical; \code{TRUE} if Newton–Raphson provided a better solution than Simulated Annealing.}
#'   \item{loglik}{Final log-likelihood value.}
#'   \item{aic}{Akaike Information Criterion.}
#'   \item{bic}{Bayesian Information Criterion.}
#'   \item{family}{The model family object used for fitting.}
#' }
#'
#' @importFrom stats logLik AIC BIC
#'
#' @export
summary.gnlmsa <- function(object, level = 0.95, test = c("Wald", "Rao", "LRT"), expected = TRUE, ...) {
  test <- match.arg(test)
  if (test != "Wald") test <- "Wald"

  family <- object$family
  npar_mu <- object$npar_mu
  beta <- object$beta
  gamma <- object$gamma
  one_parameter <- object$one_parameter

  v <- diag(vcov.gnlmsa(object, expected))

  se_beta <- sqrt(v[1:npar_mu])

  se_gamma <- sqrt(v[(npar_mu + 1):length(v)])

  ci <- confint.gnlmsa(object, level = level, test = test, expected = expected)

  ci_beta <- ci[1:npar_mu, ]
  ci_gamma <- ci[(npar_mu + 1):length(v), ]

  coefficients.mu <- data.frame(
    est = beta,
    se = se_beta,
    ci_beta
  )

  if (one_parameter) {
    coefficients.phi <- data.frame(
      est = 1,
      se = NA,
      NA, NA
    )
    rownames(coefficients.phi) <- "phi"
  } else {
    coefficients.phi <- data.frame(
      est = gamma,
      se = se_gamma,
      ci_gamma
    )
  }



  sig.level <- 1 - level

  ci_names <- paste0(c("inf", "sup"), "_", paste0(c(sig.level / 2 * 100, (1 - sig.level / 2) * 100), "%"))
  colnames(coefficients.mu) <- colnames(coefficients.phi) <- c("Estimate", "SE", ci_names)


  out <- list(
    coefficients.mu = coefficients.mu,
    coefficients.phi = coefficients.phi,
    nr_failed = object$nr_failed,
    nr_better = object$nr_better,
    loglik = as.numeric(logLik(object)),
    aic = AIC(object),
    bic = BIC(object),
    family = family,
    test = test,
    one_parameter = one_parameter
  )
  class(out) <- "summary.gnlmsa"
  out
}



#-------------------------------------------------------------------------------



#' Print Method for Summary of Generalized Non-Linear Models
#'
#' Displays a formatted summary for a fitted Generalized Non-Linear Model (GNLM) object of class `"summary.gnlmsa"`.
#' The output includes the model family, the optimization status, the estimated coefficients with their
#' standard errors and confidence intervals, and key model fit statistics.
#'
#' @param x An object of class `"summary.gnlmsa"`, typically returned by [summary.gnlmsa()].
#' @param digits Integer. The number of significant digits to use when printing.
#'   Defaults to \code{max(3L, getOption("digits") - 3L)}.
#' @param ... Further arguments passed to or from other methods (currently unused).
#'
#' @details
#' The printed summary includes:
#' \itemize{
#'   \item The model family used.
#'   \item The status of the Newton–Raphson optimization (success or fallback to Simulated Annealing).
#'   \item Estimated coefficients for the mean and dispersion components, along with standard errors
#'         and confidence intervals.
#'   \item The type of test used for confidence intervals (currently only Wald).
#'   \item Log-likelihood, Akaike Information Criterion (AIC), and Bayesian Information Criterion (BIC) values.
#' }
#'
#' @return The object \code{x} is returned invisibly.
#' @importFrom stats printCoefmat
#' @export
print.summary.gnlmsa <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nGeneralized Non-Linear Model")
  cat(paste0("\nFamily: ", x$family$family), "\n")

  if (!x$nr_failed) {
    cat("\nNewton-Raphson optimization: SUCCEEDED\nReporting Newton-Raphson estimates\n")
  } else {
    cat("\nNewton-Raphson optimization: FAILED\nReporting Simulated Annealing estimates\n")
  }

  cat("\nCoefficients for mean component:\n")
  printCoefmat(x$coefficients.mu, digits = digits)
  cat("\nCoefficients for dispersion component", if (x$one_parameter) " (fixed):\n" else ":\n", sep = "")
  printCoefmat(x$coefficients.phi, digits = digits)
  cat(paste0("\nConfidence interval: ", x$test, "\n"))
  cat("\n")
  print.default(c(loglik = x$loglik, AIC = x$aic, BIC = x$bic), digits = digits)
  invisible(x)
}



#-------------------------------------------------------------------------------



#' Residual method for \code{gnlmsa}
#'
#' Compute residuals for objects of class \code{gnlmsa}.
#'
#' @param object An object of class \code{gnlmsa}, typically returned by a call to the fitting function (e.g., \code{gnlmsa()}).
#' @param type Character. The type of residuals to compute. One of:
#'   \describe{
#'     \item{\code{"pearson"}}{Pearson residuals: \eqn{(y_i - \mu_i) / \sqrt{v(\mu_i, \phi_i)}}.}
#'     \item{\code{"response"}}{Response residuals: \eqn{y - \mu_i}.}
#'     \item{\code{"deviance"}}{Deviance residuals: \eqn{\mathrm{sign}(y - \mu)\sqrt{2\{\,\ell(y_i; y_i, \phi_i) - \ell(y_i; \mu_i, \phi_i)\}}}, where \eqn{\ell(\cdot)} is the log-likelihood.}
#'   }
#'   Defaults to \code{c("pearson", "response", "deviance")}.
#' @param ...  further arguments passed to or from other methods..
#'
#' @return A numeric vector of residuals of the specified type, one for each observation.
#'
#' @details
#' This S3 method dispatches on class \code{gnlmsa} and computes residuals using
#' the fitted mean \eqn{\mu}, dispersion \eqn{\phi}, and variance function
#' \eqn{v(\mu, \phi)} stored in the \code{object}.
#'
#' @export
residuals.gnlmsa <- function(object, type = c("pearson", "response", "deviance"), ...) {
  y <- object$y
  mu <- object$mu
  phi <- object$phi
  v <- object$family$variance(mu, phi)
  type <- match.arg(type)

  if (type == "pearson") {
    (y - mu) / sqrt(v)
  } else if (type == "response") {
    y - mu
  } else if (type == "deviance") {
    l <- object$family$loglik
    sign(y - mu) * sqrt(2 * (l(y, y, phi) - l(y, mu, phi)))
  }
}










