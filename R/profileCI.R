#' Confidence Intervals using Profile Likelihood
#'
#' Calculates confidence intervals for one or more parameters in a fitted
#' model object. A function that returns the log-likelihood must be supplied,
#' either directly via the argument `loglik` or using a [`logLikFn`] S3
#' generic.
#'
#' @param object A fitted model object. This object must have a `coef` S3
#'   method. If `faster = TRUE` then it must also have a `vcov` S3 method.
#' @param loglik A named function that returns the log-likelihood based on
#'   input parameter values and data. The first argument must be the vector of
#'   model parameters. If the likelihood is zero for any observation in the
#'   data then the function should return `-Inf.`
#'
#'   Alternatively, `loglik` does not need to be supplied if a [`logLikFn`] S3
#'   method has been created for `object`. The `profileCI` package provides
#'   `logLikFn.glm`, which is used in an example in **Examples**.
#' @param ... Further arguments to be passed to `loglik`.
#' @param parm A vector specifying the parameters for which confidence
#'   intervals are calculated, either a vector of numbers or a vector of names.
#'   The default, `which = "all"`, produces confidence intervals for all the
#'   parameters.
#' @param level The confidence level required.  A numeric scalar in (0, 1).
#' @param profile A logical scalar. If `TRUE` then confidence intervals
#'   based on a profile log-likelihood are returned.  If `FALSE` then intervals
#'   based on approximate large sample normal theory, which are symmetric about
#'   the MLE, are returned.
#' @param mult A positive numeric scalar. Controls the increment by which the
#'   parameter of interest is increased/decreased when profiling above/below
#'   its MLE. The increment is `mult * se / 100` where `se` is the estimated
#'   standard error of the estimator of the parameter. Decreasing `mult`
#'   profiles at more points but will be slower.
#' @param faster A logical scalar. If `faster = TRUE` then the profiling of the
#'   log-likelihood is in search of a lower (upper) confidence limit is
#'   started at the corresponding symmetric lower (upper) confidence limit.
#' @param epsilon Only relevant if `profile = TRUE`. A numeric vector of values
#'   that determine the accuracy of the confidence limits. `epsilon` is
#'   recycled to the length of the parameter vector `parm`.
#'
#'   * If `epsilon[i] > 0` then this value is passed as the argument `epsilon`
#'     to the [`itp::itp`] function, which estimates the parameter values for
#'     which the profile log-likelihood for parameter `i` drops to the value
#'     that defines the confidence limits, once profiling has been successful
#'     in finding an interval within which this value lies.
#'
#'    * If `epsilon[i] < 0` quadratic interpolation is used, which will tend to
#'      be faster.
#'
#'   * If `epsilon[i] = 0` then linear interpolation is used, which will be
#'     faster still.
#' @param optim_args A list of further arguments (other than `par` and `fn`) to
#'   pass to `[stats::optim]`.
#' @details The default, `epsilon = -1`, should work well enough in most
#'   circumstances, but to achieve a specific accuracy set `epsilon` to be
#'   a small positive value, for example, `epsilon = 1e-4`.
#' @seealso [`plot.profileCI`] and [`print.profileCI`].
#' @examples
#' ## From example(glm)
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3, 1, 9)
#' treatment <- gl(3, 3)
#' glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())
#' confint(glm.D93)
#' confint.default(glm.D93)
#'
#' poisson_loglik <- function(pars) {
#'   lambda <- exp(model.matrix(glm.D93) %*% pars)
#'   loglik <- stats::dpois(x = glm.D93$y, lambda = lambda, log = TRUE)
#'   return(sum(loglik))
#' }
#'
#' # Will be slower than profile.glm() because glm.fit() is fast
#' x <- profileCI(glm.D93, loglik = poisson_loglik, mult = 32)
#' x
#' plot(x, parm = 1)
#' plot(x, parm = "outcome2")
#'
#' # A logLikFn.glm S3 method is provided in profileCI so we do not need to
#' # supply loglik explicitly
#' x <- profileCI(glm.D93, mult = 32)
#' x
#'
#' x <- profileCI(glm.D93, loglik = poisson_loglik, mult = 32, faster = TRUE)
#' x
#' @export
profileCI <- function(object, loglik, ..., parm = "all", level = 0.95,
                      profile = TRUE, mult = 2, faster = FALSE, epsilon = -1,
                      optim_args = list()) {
  # Check that the model has at least 2 parameters
  cf <- coef(object)
  if (length(cf) < 2) {
    stop("The model must have more than one parameter.")
  }
  # If loglik is missing then check whether object has a logLikFn method
  # If it does then use it, otherwise throw an error
  if (missing(loglik)) {
    find_logLikFn <- function(x) {
      return(!is.null(getS3method("logLikFn", x, optional = TRUE)))
    }
    has_logLikFnMethod <- sapply(class(object), FUN = find_logLikFn)
    if (any(has_logLikFnMethod)) {
      loglik <- function(pars) {
        return(logLikFn(object, pars = pars))
      }
    } else {
      stop("")
    }
  }
  # Check and set parm
  if (is.null(names(cf))) {
    names(cf) <- paste0("par", 1:length(cf))
  }
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

  # Use the vcov method to calculate symmetric CIs
  # If profile = FALSE then we return these
  # If faster = TRUE then we pass these intervals to faster_profile_ci()
  if (!profile | faster) {
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
  # If profile = FALSE then return the symmetric intervals
  if (!profile) {
    rownames(ci_mat) <- parm
    attr(ci_mat, "interval_type") <- "symmetric"
    attr(ci_mat, "level") <- level
    class(ci_mat) <- c("profileCI", "matrix", "array")
    return(ci_mat)
  }

  # The number of parameters
  n_parm <- length(parm)
  # Force epsilon to have length equal to the number of parameters
  epsilon <- rep_len(epsilon, n_parm)
  # Extract the parameter numbers to include
  parm_numbers <- (1:length(parm_names))[which_parm]
  # An empty list in which to store the profile log-likelihood values
  # Likewise for the values of the parameters at the confidence limits
  for_plot <- list()
  lower_pars <- list()
  upper_pars <- list()
  # Calculate the estimated standard errors (to set inc below)
  ses <- sqrt(diag(vcov(object)))[parm]
  # Set up the negated log-likelihood function
  negated_loglik_fn <- function(parm, ...) {
    return(-loglik(parm, ...))
  }
  # Loop over all parameters
  # It is possible for stats::option() to throw an error, particularly if mult
  # is large. This may reflect poor starting values. in an attempt to avoid
  # this problem, we reduce (halve) mult iteratively down to the minimum value
  # min_mult. If no fit is successful then we return NAs.
  min_mult <- 1
  save_mult <- mult
  for (i in 1:n_parm) {
    success <- FALSE
    if (faster) {
      while (mult >= min_mult) {
        # Set inc based on the estimated standard errors
        inc <- mult * ses / 100
        conf_list <- faster_profile_ci(object = object,
                                       negated_loglik_fn = negated_loglik_fn,
                                       which = parm_numbers[i],
                                       which_name = parm[i],
                                       level = level, mle = coef(object),
                                       ci_sym_mat = ci_sym_mat,
                                       inc = inc[i], epsilon = epsilon[i], ...)
        if (!is.null(conf_list$optim_error)) {
          mult <- mult / 2
        } else {
          success <- TRUE
          mult <- min_mult - 1
        }
      }
    } else {
      while (mult >= min_mult) {
        # Set inc based on the estimated standard errors
        inc <- mult * ses / 100
        conf_list <- profile_ci(object = object,
                                negated_loglik_fn = negated_loglik_fn,
                                which = parm_numbers[i], level = level,
                                mle = coef(object), inc = inc[i],
                                epsilon = epsilon[i], ...)
        if (!is.null(conf_list$optim_error)) {
          mult <- mult / 2
        } else {
          success <- TRUE
          mult <- min_mult - 1
        }
      }
    }
    if (success) {
      ci_mat[i, ] <- conf_list$par_prof[c(1, 3)]
      colnames(conf_list$for_plot)[1] <- paste0(parm[i], "_values")
      for_plot[[i]] <- conf_list$for_plot
      lower_pars[[i]] <- conf_list$lower_pars
      upper_pars[[i]] <- conf_list$upper_pars
    } else {
      ci_mat[i, ] <- matrix(NA, 1, 2)
      for_plot[[i]] <- NA
      lower_pars[[i]] <- NA
      upper_pars[[i]] <- NA
    }
    # Reset mult
    mult <- save_mult
  }
  names(for_plot) <- parm

  # Format the matrix to be returned
  low <- paste0(100 * (1 - level)/ 2, "%")
  up <- paste0(100 - 100 * (1 - level)/ 2, "%")
  colnames(ci_mat) <- c(low, up)
  rownames(ci_mat) <- parm
  attr(ci_mat, "for_plot") <- for_plot
  attr(ci_mat, "crit") <- conf_list$crit
  attr(ci_mat, "interval_type") <- "profile"
  attr(ci_mat, "level") <- level
  if (any(epsilon > 0) || !faster) {
    prof_names <- paste0("when profiling for ", parm_names)
    names(lower_pars) <- prof_names[which_parm]
    names(upper_pars) <- prof_names[which_parm]
  }
  attr(ci_mat, "lower_pars") <- lower_pars
  attr(ci_mat, "upper_pars") <- upper_pars
  class(ci_mat) <- c("profileCI", "matrix", "array")
  return(ci_mat)
}
