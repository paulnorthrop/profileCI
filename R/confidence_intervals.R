#'
#'
#' @export
confint.evmiss <- function(object, parm = "all", level = 0.95, profile = FALSE,
                           mult = 2, faster = FALSE, ...) {
  # Check inputs
  parm_values <- c("mu", "sigma", "xi")
  check_values <- c("mu", "sigma", "xi", "all")
  p_message <- "c(''mu'', ''sigma'', ''xi'')"
  if (!all(is.element(parm, check_values))) {
    stop(paste("''parm'' must be ''all'', or a subset of", p_message))
  }
  if (length(parm) == 1 && parm == "all") {
    parm <- parm_values
  } else {
    parm <- parm[order(parm)]
  }
  # Logical vector indicating which parameters to include
  which_parm <- is.element(parm_values, parm)
  # Check the input confidence level
  if (level <= 0 | level >= 1) {
    stop("''level'' must be in (0, 1)")
  }

  # Calculate symmetric intervals

  z_val <- stats::qnorm(1 - (1 - level) / 2)
  mles <- coef(object)[which_parm]
  ses <- sqrt(diag(vcov(object)))[which_parm]
  sym_lower <- mles - z_val * ses
  sym_upper <- mles + z_val * ses
  ci_mat <- cbind(sym_lower, sym_upper)
  rownames(ci_mat) <- parm
  # Constrain any limits for sigma to [0, Infty]
  if (is.element("sigma", rownames(ci_mat))) {
    sigma_row <- which(rownames(ci_mat) == "sigma")
    ci_mat[sigma_row, 1] <- max(0, ci_mat[sigma_row, 1])
  }
  ci_sym_mat <- ci_mat

  # If profile log-likelihood-based intervals are required then calculate them

  if (profile) {
    # The number of parameters
    n_parm <- length(parm)
    # Extract the parameter numbers to include
    parm_numbers <- (1:length(parm_values))[which_parm]
    # Recreate the list maxima_notNA
    maxima_notNA <- list(maxima = object$maxima, notNA = object$notNA,
                         b = object$b)
    # An empty list in which to store the profile log-likelihood values
    for_plot <- list()
    # Set inc based on the estimated standard errors
    ses <- sqrt(diag(vcov(object)))
    inc <- mult * ses / 100
    # Loop over all parameters
    for (i in 1:n_parm) {
      if (faster) {
        conf_list <- faster_profile_ci(negated_loglik_fn = negated_gev_loglik,
                                       which = parm_numbers[i], level = level,
                                       mle = coef(object),
                                       ci_sym_mat = ci_sym_mat,
                                       inc = inc[i],
                                       maxima_notNA = maxima_notNA,
                                       adjust = object$adjust)
      } else {
        conf_list <- profile_ci(negated_loglik_fn = negated_gev_loglik,
                                which = parm_numbers[i], level = level,
                                mle = coef(object), inc = inc[i],
                                maxima_notNA = maxima_notNA,
                                adjust = object$adjust)
      }
      ci_mat[i, ] <- conf_list$par_prof[c(1, 3)]
      colnames(conf_list$for_plot)[1] <- paste0(parm[i], "_values")
      for_plot[[i]] <- conf_list$for_plot
    }
    names(for_plot) <- parm
  }
  # Format the matrix to be returned
  low <- paste0(100 * (1 - level)/ 2, "%")
  up <- paste0(100 - 100 * (1 - level)/ 2, "%")
  colnames(ci_mat) <- c(low, up)
  # If profile = TRUE save the profile log-likelihood values for plotting
  if (profile) {
    attr(ci_mat, "for_plot") <- for_plot
    attr(ci_mat, "crit") <- conf_list$crit
    attr(ci_mat, "level") <- level
  }
  class(ci_mat) <- c("confint_gev", "evmiss")
#  class(out) <- c("profile.profileCI", "profile")
  return(ci_mat)
}
