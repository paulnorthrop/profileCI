#' profileCI: Profiling a Log-Likelihood to Calculate Confidence Intervals
#'
#' Provides tools for profiling a user-supplied log-likelihood function to
#' calculate confidence intervals for model parameters. Speed of computation
#' can be improved by adjusting the step sizes in the profiling and/or starting
#' the profiling from limits based on the approximate large sample normal
#' distribution for the maximum likelihood estimator of a parameter. The
#' accuracy of the limits can be set by the user. A plot method visualises the
#' log-likelihood and confidence interval. Only convex log-likelihoods are
#' supported, that is, disjoint confidence intervals will not be found.
#'
#' @details The main function is [`profileCI`], which profiles the
#'   log-likelihood, with respect to one parameter at a time.
#' @docType package
#' @aliases profileCI-package
#' @seealso [`profileCI`] and [`plot.profileCI`].
#' @import stats
"_PACKAGE"
