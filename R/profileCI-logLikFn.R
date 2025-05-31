#' Calculate Log-Likelihood
#'
#' This is a generic function for calculating a log-likelihood for an object
#' for input parameter values.
#'
#' @param object A fitted model object.
#' @param pars A numeric vector of parameters of the model.
#' @param ... Further arguments.
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
