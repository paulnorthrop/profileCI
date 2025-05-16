#' Confidence Intervals using Profile Likelihood
#'
#' Calculates confidence intervals for one or more parameters for a fitted
#' model object. A function that returns the log-likelihood must also be
#' supplied.
#' Also check `chandwich::adjust_loglik` for `alg_deriv` and `alg_hess`.
#' Allow `loglik` as an attribute to object.
#' If only `loglik` provided then find MLEs etc
#' Use lower and upper, if present, to contrain symmetric intervals
#' Optional function to calculate initial estimates, conditional on 1 parameter
#'
#' @param object A fitted model object. This object must have a `coef` S3
#'   method. If `faster = TRUE` then it must also have a `vcov` S3 method.
#' @param loglik A named function that returns the log-likelihood based on
#'   input parameter values and data. The first argument must be the vector of
#'   model parameters. If the likelihood is zero for any observation in the
#'   data then the function should return `-Inf.`
#' @param ... Further arguments to be passed to `loglik`.
#' @param parm A character vector specifying the parameters for which
#'   confidence intervals are calculated, either a vector of numbers or a vector
#'   of names. The default, `which = "all"`, produces confidence intervals for
#'   all the parameters.
#' @param level The confidence level required.  A numeric scalar in (0, 1).
#' @param mult A positive numeric scalar. Controls the increment by which the
#'   parameter of interest is increased/decreased when profiling above/below
#'   its MLE. The increment is `mult * se / 100` where `se` is the estimated
#'   standard error of the estimator of the parameter. Decreasing `mult`
#'   profiles at more points but will be slower.
#' @param faster A logical scalar. If `faster = TRUE` then the profiling of the
#'   log-likelihood is in search of a lower (upper) confidence limit is
#'   started at the corresponding symmetric lower (upper) confidence limit.
#' @param epsilon Only relevant if `profile = TRUE`. A numeric vector of values
#'   passed as the argument `epsilon` to [`itp::itp`] to set the desired
#'   accuracy of the confidence limits. `epsilon` is recycled to the length of
#'   the parameter vector `parm`. If `epsilon[i]` is positive then the
#'   [`itp::itp`] function is used to estimate the parameter values for which
#'   the profile log-likelihood for parameter `i` drops to the value that
#'   defines the confidence limits, once profiling has been successful in
#'   finding an interval within which this value lies.
#'   Otherwise, (if `epsilon <= 0`) linear interpolation is used.
#' @param optim_args A list of further arguments (other than `par` and `fn`) to
#'   pass to `[stats::optim]`.
#' @details [lm loglik](https://stats.stackexchange.com/questions/73196/recalculate-log-likelihood-from-a-simple-r-lm-model)
#'   **(Use lm as an example: prof = sym)**
#'   **(Also include an example where this isn't the case: GEV?)**
#'   **(glm's confint method uses profiling)**
#'   **(Could I use ... for both arguments to loglik and optim? See chandwich::adjust_loglik())**
#'   **(It may be easier not to do this)**
#'   **(If faster = TRUE, use optim arguments lower and upper to constrain the symmetric limits)**
#' @examples
#' ## From example(glm)
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3, 1, 9); treatment <- gl(3, 3)
#' glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())
#' confint(glm.D93)
#' confint.default(glm.D93)
#'
#' poisson_loglik <- function(par) {
#'   lambda <- exp(model.matrix(glm.D93) %*% par)
#'   loglik <- stats::dpois(x = counts, lambda = lambda, log = TRUE)
#'   return(sum(loglik))
#' }
#' # Will be slower than profile.glm() because glm.fit() is fast
#' x <- profileCI(glm.D93, loglik = poisson_loglik, mult = 32)
#' x
#' plot(x, parm = 1)
#' plot(x, parm = "outcome2")
#'
#' # Can I speed things up by starting at the symmetric limits?
#' # Do I need to estimate other parameters conditional on the parameter of interest?
#' # Base this on the vcov matrix
#' @export
profileCI <- function(object, loglik, ..., parm = "all", level = 0.95,
                      mult = 2, faster = FALSE, epsilon = 1e-4,
                      optim_args = list()) {

  # Check and set parm
  cf <- coef(object)
  parm_names <- names(cf)
  if (is.character(parm)) {
    if (!all(is.element(parm, c(parm_names, "all")))) {
      p_message <- paste0("''", parm_names, "''", collapse = ",")
      stop("Character, ''parm'' must be ''all'', or a subset of ", p_message)
    }
  } else {
    if (!all(is.element(parm, 1:length(cf)))) {
      p_message <- paste0(1:length(cf), collapse = ",")
      stop("Numeric, ''parm'' must be ''all'', or a subset of ", p_message)
    }
  }
  if (length(parm) == 1 && parm == "all") {
    parm <- parm_names
  } else if (is.numeric(parm)) {
    parm <- parm_names[parm]
  }
  # Logical vector indicating which parameters to include
  which_parm <- is.element(parm_names, parm)
  # Check the input confidence level
  if (level <= 0 | level >= 1) {
    stop("''level'' must be in (0, 1)")
  }

  # If faster = TRUE then we need the vcov method to calculate symmetric CIs
  if (faster) {
    z_val <- stats::qnorm(1 - (1 - level) / 2)
    mles <- coef(object)[parm]
    ses <- sqrt(diag(vcov(object)))[parm]
    sym_lower <- mles - z_val * ses
    sym_upper <- mles + z_val * ses
    ci_mat <- cbind(sym_lower, sym_upper)
    rownames(ci_mat) <- parm
    low <- paste0(100 * (1 - level)/ 2, "%")
    up <- paste0(100 - 100 * (1 - level)/ 2, "%")
    colnames(ci_mat) <- c(low, up)
    ci_sym_mat <- ci_mat
  } else {
    ci_mat <- matrix(NA, ncol = 2, nrow = length(parm))
  }

  # The number of parameters
  n_parm <- length(parm)
  # Extract the parameter numbers to include
  parm_numbers <- (1:length(parm_names))[which_parm]
  # An empty list in which to store the profile log-likelihood values
  for_plot <- list()
  # Set inc based on the estimated standard errors
  ses <- sqrt(diag(vcov(object)))
  inc <- mult * ses / 100
  # Set up the negated log-likelihood function
  negated_loglik_fn <- function(parm, ...) {
    return(-loglik(parm, ...))
  }
  # Loop over all parameters
  for (i in 1:n_parm) {
    if (faster) {
      conf_list <- faster_profile_ci(negated_loglik_fn = negated_loglik_fn,
                                     which = parm_numbers[i], level = level,
                                     mle = coef(object),
                                     ci_sym_mat = ci_sym_mat,
                                     inc = inc[i], epsilon = epsilon, ...)
    } else {
      conf_list <- profile_ci(negated_loglik_fn = negated_loglik_fn,
                              which = parm_numbers[i], level = level,
                              mle = coef(object), inc = inc[i],
                              epsilon = epsilon, ...)
    }
    ci_mat[i, ] <- conf_list$par_prof[c(1, 3)]
    colnames(conf_list$for_plot)[1] <- paste0(parm[i], "_values")
    for_plot[[i]] <- conf_list$for_plot
  }
  names(for_plot) <- parm

  # Format the matrix to be returned
  low <- paste0(100 * (1 - level)/ 2, "%")
  up <- paste0(100 - 100 * (1 - level)/ 2, "%")
  colnames(ci_mat) <- c(low, up)
  rownames(ci_mat) <- parm_names
  attr(ci_mat, "for_plot") <- for_plot
  attr(ci_mat, "crit") <- conf_list$crit
  attr(ci_mat, "level") <- level
  class(ci_mat) <- c("profileCI", "matrix", "array")
  return(ci_mat)
}
