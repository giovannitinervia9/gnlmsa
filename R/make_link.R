


#-------------------------------------------------------------------------------



#' Create a link function object for generalized non-linear models
#'
#' This function creates a link function object similar to `stats::make.link()` but extends it by
#' including the second derivative of the mean with respect to eta (`mu2.eta2`).
#' This is useful for generalized non-linear models where higher-order derivatives are required.
#'
#' @param link Character string specifying the link function. Options include:
#'   \itemize{
#'     \item \code{"logit"}: log(p/(1-p))
#'     \item \code{"probit"}: inverse of the standard normal CDF
#'     \item \code{"cauchit"}: inverse of the Cauchy CDF
#'     \item \code{"cloglog"}: log(-log(1-p))
#'     \item \code{"identity"}: no transformation
#'     \item \code{"log"}: natural logarithm
#'     \item \code{"sqrt"}: square root
#'     \item \code{"1/mu^2"}: inverse square
#'     \item \code{"inverse"}: reciprocal
#'   }
#'
#' @return A list with S3 class \code{"link-gnlm"} containing the following components:
#'   \itemize{
#'     \item \code{linkfun}: Function that transforms the mean (mu) to eta
#'     \item \code{linkinv}: Function that transforms eta to the mean (mu)
#'     \item \code{mu.eta}: Function returning the derivative of mu with respect to eta (\eqn{d\mu/d\eta})
#'     \item \code{mu2.eta2}: Function returning the second derivative of mu with respect to eta (\eqn{d^2\mu/d^2\eta^2})
#'     \item \code{valideta}: Function that checks if the eta values are valid
#'     \item \code{name}: Character string naming the link function
#'   }
#'
#' @note This function extends \code{stats::make.link()} by including the \code{mu2.eta2} component,
#' which calculates the second derivative of the mean with respect to eta.
#' While in standard GLMs eta is typically a linear predictor, in non-linear models eta can be
#' a non-linear function of the parameters, which is why these functions are particularly useful
#' for generalized non-linear models.
#'
#' @examples
#' # Create a logit link object
#' logit_link <- make_link("logit")
#'
#' # Apply link functions to a probability value
#' mu <- 0.75
#' eta <- logit_link$linkfun(mu)  # Transform to eta
#' mu_back <- logit_link$linkinv(eta)  # Transform back to mean
#'
#' # Calculate derivatives at eta
#' first_deriv <- logit_link$mu.eta(eta)
#' second_deriv <- logit_link$mu2.eta2(eta)
#'
#' @export
#' @importFrom stats dcauchy dnorm pcauchy pnorm qcauchy qnorm
make_link <- function (link) {
  switch(link,

         logit = {
           linkfun <- function(mu) .Call(C_logit_link, mu)
           linkinv <- function(eta) .Call(C_logit_linkinv, eta)
           mu.eta <- function(eta) .Call(C_logit_mu_eta, eta)
           mu2.eta2 <- function(eta) 2*exp(2*eta)/(1 + exp(eta))^3 - exp(eta)/(1 + exp(eta))^2
           valideta <- function(eta) TRUE
         },

         probit = {
           linkfun <- function(mu) qnorm(mu)
           linkinv <- function(eta) {
             thresh <- -qnorm(.Machine$double.eps)
             eta <- pmin(pmax(eta, -thresh), thresh)
             pnorm(eta)
           }
           mu.eta <- function(eta) pmax(dnorm(eta), .Machine$double.eps)
           mu2.eta2 <- function(eta) pmax(eta*dnorm(eta), .Machine$double.eps)
           valideta <- function(eta) TRUE
         },

         cauchit = {
           linkfun <- function(mu) qcauchy(mu)
           linkinv <- function(eta) {
             thresh <- -qcauchy(.Machine$double.eps)
             eta <- pmin(pmax(eta, -thresh), thresh)
             pcauchy(eta)
           }
           mu.eta <- function(eta) pmax(dcauchy(eta), .Machine$double.eps)
           mu2.eta2 <- function(eta) pmax(-2*eta/(pi*(1 + eta^2)^2), .Machine$double.eps)
           valideta <- function(eta) TRUE
         },

         cloglog = {
           linkfun <- function(mu) log(-log(1 - mu))
           linkinv <- function(eta) pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
           mu.eta <- function(eta) {
             eta <- pmin(eta, 700)
             pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
           }
           mu2.eta2 <- function(eta){
             eta <- pmin(eta, 700)
             pmax(-exp(eta - exp(eta))*(exp(eta) - 1) , .Machine$double.eps)
           }
           valideta <- function(eta) TRUE
         },

         identity = {
           linkfun <- function(mu) mu
           linkinv <- function(eta) eta
           mu.eta <- function(eta) rep.int(1, length(eta))
           mu2.eta2 <- function(eta) rep.int(0, length(eta))
           valideta <- function(eta) TRUE
         },

         log = {
           linkfun <- function(mu) log(mu)
           linkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)
           mu.eta <- function(eta) pmax(exp(eta), .Machine$double.eps)
           mu2.eta2 <- function(eta) pmax(exp(eta), .Machine$double.eps)
           valideta <- function(eta) TRUE
         },

         sqrt = {
           linkfun <- function(mu) sqrt(mu)
           linkinv <- function(eta) eta^2
           mu.eta <- function(eta) 2 * eta
           mu2.eta2 <- function(eta) 2
           valideta <- function(eta) all(is.finite(eta)) && all(eta > 0)
         },

         `1/mu^2` = {
           linkfun <- function(mu) 1/mu^2
           linkinv <- function(eta) 1/sqrt(eta)
           mu.eta <- function(eta) -1/(2 * eta^1.5)
           mu2.eta2 <- function(eta) 3/(4 * eta^2.5)
           valideta <- function(eta) all(is.finite(eta)) && all(eta > 0)
         },

         inverse = {
           linkfun <- function(mu) 1/mu
           linkinv <- function(eta) 1/eta
           mu.eta <- function(eta) -1/(eta^2)
           mu2.eta2 <- function(eta) 2/eta^3
           valideta <- function(eta) all(is.finite(eta)) && all(eta != 0)
         },
         stop(gettextf("%s link not recognised", sQuote(link)),
              domain = NA))
  environment(linkfun) <- environment(linkinv) <- environment(mu.eta) <- environment(valideta) <- asNamespace("stats")
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, mu2.eta2 = mu2.eta2,
                 valideta = valideta, name = link), class = "link-gnlm")
}




