#' Calculate Log-Likelihood
#'
#' This is a generic function for calculating a log-likelihood for an object
#' at input parameter values.
#'
#' @param object A fitted model object.
#' @param pars A numeric vector of parameters of the model.
#' @param ... Further arguments.
#' @details This generic function has been created to enable a function that
#'   calculates the log-likelihood for a parametric model at a given set of
#'   parameter values in `pars` to be available to the function [`profileCI`].
#' @return A numeric scalar. The value of the log-likelihood function for the
#'   fitted model object `object` for parameter values `pars`.
#'
#'   The `logLikFn.glm` generic is specifically for the Poisson log-linear GLM
#'   case.
#' @name logLikFn
NULL

#' @rdname logLikFn
#' @export
logLikFn <- function(object, pars, ...) {
  UseMethod("logLikFn")
}

#' @rdname logLikFn
#' @export
logLikFn.glm <- function(object, pars, ...) {
  lambda <- exp(model.matrix(object) %*% pars)
  loglik <- stats::dpois(x = object$y, lambda = lambda, log = TRUE)
  return(sum(loglik))
}
