#' Methods for objects of class `"profileCI"`
#'
#' Methods for objects of class `"profileCI"` returned from [`profileCI`].
#' @param x An object inheriting from class `profileCI"`, a result of a call
#'   to [`profileCI`].
#' @param ... Further arguments. For `print."profileCI"` to pass arguments to
#'   [`print`]. For `plot.profileCI` to pass graphical parameters to
#'   [`graphics::plot`] to create the initial plot of the profile
#'   log-likelihood.
#' @details `print.profileCI`. A numeric matrix with 2 columns giving the
#'   lower and upper confidence limits for the parameters specified by the
#'   argument `parm` in [`profileCI`]. These columns are labelled as
#'   `(1-level)/2` and `1-(1-level)/2`, expressed as a percentage, by default
#'   `2.5%` and `97.5%`.
#'
#'   `plot.profileCI`. A plot is produced of the profile log-likelihood for
#'   the parameter chosen by `parm`. Only the parameter values used to profile
#'   the log-likelihood in the call to [`profileCI`] are included, so
#'   if `faster = TRUE` was used then the plot will not be of a smooth curve
#'   but will be triangular in the middle.
#' @return `print.profileCI`: the argument `x` is returned, invisibly.
#'
#'   `plot.profileCI`: a numeric vector containing the confidence interval
#'   for the parameter chosen for the plot.
#' @section Examples: See [`profileCI`].
#' @seealso [`profileCI`].
#' @name profileCI_methods
NULL
## NULL

# ============================= print.profileCI ============================= #

#' Print method for objects of class `"profileCI"`
#'
#' @rdname profileCI_methods
#' @export
print.profileCI <- function(x, ...) {
  print(x[1:nrow(x), , drop = FALSE], ...)
  return(invisible(x))
}

# ============================== plot.profileCI ============================= #

#' Plot method for objects of class `"profileCI"`
#'
#' @param parm A numeric or character scalar specifying the parameter for which
#'   a profile log-likelihood is plotted. Must be a single component consistent
#'   with the argument `parm` to [`profileCI`].
#' @param add A logical scalar. If `add = TRUE` then the plot is annotated with
#'   a horizontal line indicating the critical value for the profile
#'   log-likelihood used to calculate the confidence limits, vertical lines
#'   indicating the values of these limits and a legend stating the
#'   confidence interval.
#' @param digits An integer. Passed to [`signif`] to round the confidence
#'   limits in the legend, if `add = TRUE`.
#' @rdname profileCI_methods
#' @export
plot.profileCI <- function(x, parm = 1:nrow(x), add = TRUE, digits = 2, ...) {
  # If symmetric intervals were produced then do not create a plot
  if (attr(x, "interval_type") == "symmetric") {
    stop("There is no plot for 'profile = FALSE' in profileCI().")
  }
  # For which parameter is the plot required?
  # Work with the parameter number rather than the parameter name
  parm_names <- rownames(x)
  n_pars <- length(parm_names)
  if (is.numeric(parm)) {
    parm <- parm[1]
    if (!any(is.element(parm, 1:n_pars))) {
      pars_message <- paste0("{", paste0(1:n_pars, collapse = ","), "}")
      stop("''If parm is numeric it must be in ''", pars_message)
    }
  } else {
    if (any(!is.element(parm, parm_names))) {
      stop(parm, " was not included in ''parm'' in the call to profileCI().")
    }
    parm <- which(is.element(parm_names, parm))
  }
  to_plot <- attr(x, "for_plot")[[parm]]
  my_xlab <- parm_names[parm]
  limits <- x[parm, ]
  crit <- attr(x, "crit")
  plot_fn <- function(x, ..., xlab = my_xlab, ylab = "profile log-likelihood",
                      lwd = 2, type = "l") {
    graphics::plot(x, ..., xlab = xlab, ylab = ylab, lwd = lwd, type = type)
  }
  plot_fn(to_plot, ...)
  if (add) {
    user_args <- list(...)
    graphics::abline(h = crit, lty = 2)
    graphics::abline(v = limits, lty = 2)
    level <- attr(x, "level") * 100
    rlimits <- signif(limits, digits)
    legend_text <- paste0(level, "% CI: (", rlimits[1], ",", rlimits[2], ")")
    graphics::legend("bottom", legend = legend_text)
  }
  return(invisible(limits))
}
