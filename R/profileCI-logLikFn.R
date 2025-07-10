#' Calculate Log-Likelihood
#'
#' This is a generic function for calculating a log-likelihood for an object
#' at input parameter values.
#'
#' @param object A fitted model object.
#' @param pars A numeric vector of parameters of the model.
#' @param ... Further arguments. None are used in the `logLikFn.glm` and
#'   `logLikFn.nls` generics.
#' @details This generic function has been created to enable a function that
#'   calculates the log-likelihood for a parametric model at a given set of
#'   parameter values in `pars` to be available to the function [`profileCI`].
#' @return A numeric scalar. The value of the log-likelihood function for the
#'   fitted model object `object` for parameter values `pars`.
#'
#'   The `logLikFn.glm` generic is specifically for the unweighted Poisson
#'   log-linear GLM case.
#'
#'   The `logLikFn.nls` generic is more general and should work for all
#'   model objects returned by [`stats::nls`].
#' @seealso [`profileCI`], [`stats::glm`], [`stats::nls`].
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

#' @rdname logLikFn
#' @export
logLikFn.nls <- function(object, pars, ...) {
  # Extract the original call and make it a list
  call_list <- as.list(object$call)
  # Use the argument start to set the parameter values
  call_list$start <- pars
  names(call_list$start) <- names(coef(object))
  # Extract the control argument to nls()
  nls_control <- call_list$control
  # Set tol = Inf so that the algorithm stops at the input parameter values
  nls_control$tol <- Inf
  # Insert into nls_control
  call_list$control <- nls_control
  # Call stats:nls
  val <- do.call(stats::nls, call_list[-1], envir = environment(object$m$conv))
  # Calculate the log-likelihood, incorporating the weights
  res <- val$m$resid()
  n <- length(res)
  w <- object$weights %||% rep_len(1, n)
  loglik <- -n * (log(2 * pi) + 1 - log(n) + log(sum(w * res ^ 2))) / 2
  return(sum(loglik))
}
