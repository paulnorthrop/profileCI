#' Internal profileCI functions
#'
#' Internal profileCI functions
#' @details
#' These functions are not intended to be called by the user.
#' @name profileCI-internal
#' @keywords internal
NULL

# ======================== General profiling function ======================= #

#' @keywords internal
#' @rdname profileCI-internal
profile_ci <- function(negated_loglik_fn, which = 1, level, mle, inc, epsilon,
                       ...) {

  # The number of parameters
  n_pars <- length(mle)

  # The -log-likelihood to profile over parameters in par other than par[which]
  profiling_fn <- function(par, par_which, ...) {
    # Place par_which in position which and the other parameters around it
    parameters <- rep_len(NA, length(par) + 1)
    parameters[which] <- par_which
    parameters[-which] <- par
    # Return the negated log-likelihood
    return(negated_loglik_fn(parameters, ...))
  }

  # The maximised log-likelihood and the MLE for the parameter of interest
  max_loglik <- -negated_loglik_fn(mle, ...)
  mle_which <- mle[which]
  # The horizontal line that determines the confidence limits
  conf_line <- max_loglik - 0.5 * stats::qchisq(level, 1)
  # Vectors to store values of the parameters (x1 and x2) and the values of
  # the profile log-likelihood (v1 and v2)
  v1 <- v2 <- x1 <- x2 <- NULL
  x2[1] <- x1[1] <- mle[which]
  v2[1] <- v1[1] <- max_loglik

  # Starting from the MLE, we search upwards and downwards until we pass the
  # cutoff for the 100level% confidence interval

  ### Upper tail ...
  par_which <- mle_which
  my_val <- max_loglik
  ii <- 1
  sol <- mle[-which]
  while (my_val > conf_line){
    par_which <- par_which + inc
    opt <- stats::optim(sol, profiling_fn, par_which = par_which, ...)
    sol <- opt$par
    ii <- ii + 1
    x2[ii] <- par_which
    v2[ii] <- -opt$value
    my_val <- v2[ii]
  }
  # Save sol for possible use later by itp()
  sol_upper <- sol
  # Also save the values of the parameters after the profile log-likelihood
  # has dropped below conf_line. These may be useful in providing initial
  # values for profiling with respect to a return level
  upper_pars <- numeric(n_pars)
  upper_pars[which] <- par_which
  upper_pars[-which] <- sol_upper
  names(upper_pars) <- names(mle)

  ### Lower tail ...
  par_which <- mle_which
  my_val <- max_loglik
  ii <- 1
  sol <- mle[-which]
  while (my_val > conf_line){
    par_which <- par_which - inc
    opt <- stats::optim(sol, profiling_fn, par_which = par_which, ...)
    sol <- opt$par
    ii <- ii + 1
    x1[ii] <- par_which
    v1[ii] <- -opt$value
    my_val <- v1[ii]
  }
  # Save sol for possible use later by itp()
  sol_lower <- sol
  # Also save the values of the parameters after the profile log-likelihood
  # has dropped below conf_line. These may be useful in providing initial
  # values for profiling with respect to a return level
  lower_pars <- numeric(n_pars)
  lower_pars[which] <- par_which
  lower_pars[-which] <- sol_lower
  names(lower_pars) <- names(mle)

  # Find the limits of the confidence interval
  prof_lik <- c(rev(v1), v2)
  par_values <- c(rev(x1), x2)

  # Find where the curve crosses conf_line
  temp <- diff(prof_lik - conf_line > 0)
  # Find the upper limit of the confidence interval
  loc <- which(temp == -1)
  x1up <- par_values[loc]
  x2up <- par_values[loc + 1]
  y1up <- prof_lik[loc]
  y2up <- prof_lik[loc + 1]
  # Find the lower limit of the confidence interval
  loc <- which(temp == 1)
  x1low <- par_values[loc]
  x2low <- par_values[loc + 1]
  y1low <- prof_lik[loc]
  y2low <- prof_lik[loc + 1]

  # If epsilon > 0 then use itp::itp()
  if (epsilon > 0) {
    itp_function <- function(x, par, ...) {
      opt <- stats::optim(par = par, fn = profiling_fn, par_which = x, ...)
      val <- -opt$value - conf_line
      attr(val, "parameters") <- opt$par
      return(val)
    }
    # Find the upper limit of the confidence interval
    upper <- itp::itp(f = itp_function, interval = c(x1up, x2up),
                      par = sol_upper,
                      f.a = y1up - conf_line, f.b = y2up - conf_line,
                      epsilon = epsilon, ...)
    up_lim <- upper$root
    # Find the lower limit of the confidence interval
    lower <- itp::itp(f = itp_function, interval = c(x1low, x2low),
                      par = sol_lower,
                      f.a = y1low - conf_line, f.b = y2low - conf_line,
                      epsilon = epsilon, ...)
    low_lim <- lower$root
    # Add the approximate solutions to par_values and prof_lik
    # Note that only the extreme ends of these vectors have prof_lik values
    # below conf_line
    n <- length(par_values)
    par_values <- c(par_values[1], low_lim, par_values[2:(n-1)], up_lim,
                    par_values[n])
    prof_lik <- c(prof_lik[1], lower$f.root + conf_line, prof_lik[2:(n-1)],
                  upper$f.root + conf_line, prof_lik[n])
    # Replace the values in lower_pars and upper_pars with those that apply to
    # the solutions found by itp::itp()
    lower_pars <- numeric(n_pars)
    lower_pars[which] <- lower$root
    lower_pars[-which] <- attr(lower$f.root, "parameters")
    names(lower_pars) <- names(mle)
    upper_pars <- numeric(n_pars)
    upper_pars[which] <- upper$root
    upper_pars[-which] <- attr(upper$f.root, "parameters")
    names(upper_pars) <- names(mle)
  } else {
    up_lim <- x1up + (conf_line - y1up) * (x2up - x1up) / (y2up - y1up)
    low_lim <- x1low + (conf_line - y1low) * (x2low - x1low) / (y2low - y1low)
    lower_pars <- NULL
    upper_pars <- NULL
  }

  par_prof <- c(lower = low_lim, mle_which, upper = up_lim)
  return(list(par_prof = par_prof, crit = conf_line,
              for_plot = cbind(par_values = par_values,
                               prof_loglik = prof_lik),
              lower_pars = lower_pars, upper_pars = upper_pars))
}

# ======================== Faster GEV profiling function ==================== #

#' @keywords internal
#' @rdname profileCI-internal
faster_profile_ci <- function(negated_loglik_fn, which = 1, which_name, level,
                              mle, ci_sym_mat, inc, epsilon, ...) {

  # The number of parameters
  n_pars <- length(mle)

  # The -log-likelihood to profile over parameters in par other than par[which]
  profiling_fn <- function(par, par_which, ...) {
    # Place par_which in position which and the other parameters around it
    parameters <- rep_len(NA, length(par) + 1)
    parameters[which] <- par_which
    parameters[-which] <- par
    # Return the negated log-likelihood
    return(negated_loglik_fn(parameters, ...))
  }

  # The maximised log-likelihood and the MLE for the parameter of interest
  max_loglik <- -negated_loglik_fn(mle, ...)
  mle_which <- mle[which]
  # The horizontal line that determines the confidence limits
  conf_line <- max_loglik - 0.5 * stats::qchisq(level, 1)
  # Vectors to store values of the parameters (x1 and x2) and the values of
  # the profile log-likelihood (v1 and v2)
  v1 <- v2 <- x1 <- x2 <- NULL
  x2[1] <- x1[1] <- mle[which]
  v2[1] <- v1[1] <- max_loglik

  ### Upper tail ...

  # We start from the upper limit of the symmetric confidence interval, using
  # ? to set initial estimates of the parameters other than the parameter which

  par_which <- ci_sym_mat[which_name, 2]
  init <- mle[-which]

  # Calculate the profile log-likelihood at the initial values
  # If this is greater than conf_line then we search upwards
  # If this is smaller than conf_line then we search downwards

  ii <- 2
  opt <- stats::optim(init, profiling_fn, par_which = par_which, ...)
  my_val <- -opt$value
  x2[ii] <- par_which
  v2[ii] <- my_val
  if (my_val < conf_line) {
    searched_upwards <- FALSE
    while_condition <- function(my_val) {
      return(my_val < conf_line)
    }
    delta <- -inc
  } else {
    searched_upwards <- TRUE
    while_condition <- function(my_val) {
      return(my_val > conf_line)
    }
    delta <- inc
  }
  sol <- opt$par
  while (while_condition(my_val)){
    par_which <- par_which + delta
    opt <- stats::optim(sol, profiling_fn, par_which = par_which, ...)
    sol <- opt$par
    ii <- ii + 1
    x2[ii] <- par_which
    v2[ii] <- -opt$value
    my_val <- v2[ii]
  }
  # If we searched downwards then reorder the results
  if (!searched_upwards) {
    x2[-1] <- rev(x2[-1])
    v2[-1] <- rev(v2[-1])
  }
  # Save sol for possible use later by itp()
  sol_upper <- sol

  ### Lower tail ...

  # We start from the upper limit of the symmetric confidence interval, using
  # ? to set initial estimates of the parameters other than the parameter which

  par_which <- ci_sym_mat[which_name, 1]
  init <- mle[-which]

  # Calculate the profile log-likelihood at the initial values
  # If this is greater than conf_line then we search downwards
  # If this is smaller than conf_line then we search upwards

  ii <- 2
  opt <- stats::optim(init, profiling_fn, par_which = par_which, ...)
  my_val <- -opt$value
  x1[ii] <- par_which
  v1[ii] <- my_val
  if (my_val < conf_line) {
    searched_downwards <- FALSE
    while_condition <- function(my_val) {
      return(my_val < conf_line)
    }
    delta <- -inc
  } else {
    searched_downwards <- TRUE
    while_condition <- function(my_val) {
      return(my_val > conf_line)
    }
    delta <- inc
  }
  sol <- opt$par
  while (while_condition(my_val)){
    par_which <- par_which - delta
    opt <- stats::optim(sol, profiling_fn, par_which = par_which, ...)
    sol <- opt$par
    ii <- ii + 1
    x1[ii] <- par_which
    v1[ii] <- -opt$value
    my_val <- v1[ii]
  }
  # If we searched upwards then reorder the results
  if (!searched_downwards) {
    x1[-1] <- rev(x1[-1])
    v1[-1] <- rev(v1[-1])
  }
  # Save sol for possible use later by itp()
  sol_lower <- sol

  # Find the limits of the confidence interval
  prof_lik <- c(rev(v1), v2)
  par_values <- c(rev(x1), x2)

  # Find where the curve crosses conf_line
  temp <- diff(prof_lik - conf_line > 0)
  # Find the upper limit of the confidence interval
  loc_upper <- which(temp == -1)
  x1up <- par_values[loc_upper]
  x2up <- par_values[loc_upper + 1]
  y1up <- prof_lik[loc_upper]
  y2up <- prof_lik[loc_upper + 1]
  # Find the lower limit of the confidence interval
  loc_lower <- which(temp == 1)
  x1low <- par_values[loc_lower]
  x2low <- par_values[loc_lower + 1]
  y1low <- prof_lik[loc_lower]
  y2low <- prof_lik[loc_lower + 1]

  # If epsilon > 0 then use itp::itp()
  if (epsilon > 0) {
    itp_function <- function(x, par, ...) {
      opt <- stats::optim(par = par, profiling_fn, par_which = x, ...)
      val <- -opt$value - conf_line
      attr(val, "parameters") <- opt$par
      return(val)
    }
    # Find the upper limit of the confidence interval
    upper <- itp::itp(itp_function, interval = c(x1up, x2up),
                      par = sol_upper,
                      f.a = y1up - conf_line, f.b = y2up - conf_line,
                      epsilon = epsilon, ...)
    up_lim <- upper$root
    # Find the lower limit of the confidence interval
    lower <- itp::itp(itp_function, interval = c(x1low, x2low),
                      par = sol_lower,
                      f.a = y1low - conf_line, f.b = y2low - conf_line,
                      epsilon = epsilon, ...)
    low_lim <- lower$root
    # Add the approximate solutions to par_values and prof_lik
    # Note that only the extreme ends of these vectors have prof_lik values
    # below conf_line
    n <- length(par_values)
    par_values <- c(par_values[1:loc_lower], low_lim,
                    par_values[(loc_lower + 1):loc_upper], up_lim,
                    par_values[(loc_upper + 1):n])
    prof_lik <- c(prof_lik[1:loc_lower], lower$f.root + conf_line,
                  prof_lik[(loc_lower + 1):loc_upper], upper$f.root + conf_line,
                  prof_lik[(loc_upper + 1):n])
    # Replace the values in lower_pars and upper_pars with those that apply to
    # the solutions found by itp::itp()
    lower_pars <- numeric(n_pars)
    lower_pars[which] <- lower$root
    lower_pars[-which] <- attr(lower$f.root, "parameters")
    names(lower_pars) <- names(mle)
    upper_pars <- numeric(n_pars)
    upper_pars[which] <- upper$root
    upper_pars[-which] <- attr(upper$f.root, "parameters")
    names(upper_pars) <- names(mle)
  } else {
    up_lim <- x1up + (conf_line - y1up) * (x2up - x1up) / (y2up - y1up)
    low_lim <- x1low + (conf_line - y1low) * (x2low - x1low) / (y2low - y1low)
    lower_pars <- NULL
    upper_pars <- NULL
  }

  par_prof <- c(lower = low_lim, mle_which, upper = up_lim)
  return(list(par_prof = par_prof, crit = conf_line,
              for_plot = cbind(par_values = par_values,
                               prof_loglik = prof_lik),
              lower_pars = lower_pars, upper_pars = upper_pars))
}

# ================== Initial estimates for GEV profiling ==================== #

#' @keywords internal
#' @rdname profileCI-internal
gev_profile_init <- function(data, mu, sigma, xi) {
  # This function provides initial estimates of two of the GEV parameters
  # (mu, sigma, xi) given a user-supplied value for the other one.

  # This may be helpful when profiling a GEV log-likelihood to find the
  # limits of a confidence interval for one of the parameters.
  # Suppose that a Wald-type symmetric confidence interval has already been
  # calculated for the parameter of interest. The limits of this interval could
  # be used to help start the search for the limits of the interval based on
  # the corresponding profile log-likelihood.

  # The estimates calculated below are based on the LRSE(EV) approach described
  # on page 99 of the first edition of the book Reiss and Thomas (1997)
  # Statistical Analysis of Extreme Values, Birkhauser, Basel, Switzerland.
  # The following uses the notation in this book.

  # We set q0 and q2 to correspond to the minimjm and maximum of the observed
  # values in the input vector data. This enables us to ensure that the
  # resulting estimates respect the constraints on the GEV parameter space
  # given the observed values in data.

  # Which of the GEV parameters has been given? There should be exactly one
  mu_given <- !missing(mu)
  sigma_given <- !missing(sigma)
  xi_given <- !missing(xi)
  number_given <- sum(mu_given, sigma_given, xi_given)
  if (number_given != 1) {
    stop("Exactly one of mu, sigma, xi must be provided.")
  }

  # Calculate quantities that will be used throughout
  n <- length(data)
  # Values of q0, q1 and q2
  q0 <- 1 / (n + 1)
  q2 <- n / (n + 1)
  a <- sqrt(log(q2) / log(q0))
  q1 <- q0 ^ a
  # The corresponding order statisticss
  n0 <- 1
  n1 <- round((n + 1) * q1)
  n2 <- n
  ns <- c(n0, n1, n2)
  qs <- c(q0, q1, q2)
  #  Extract these order statistics from the data
  x <- sort(data)[ns]

  # If xi is given then we set mu and sigma so that the fitted quantiles
  # for q0 and q2 equal the corresponding order statistics in x

  # If either mu or sigma is given then first we estimate xi. Then we use
  # either the sample minimum or maximum, depending on the sign of the estimate
  # of xi, to estimate the remaining unknown parameter.
  # If this sign is negative then we use the sample maximum.
  # If this sign is positive then we use the sample minimum.
  # This ensures that none of the data are outside the GEV support.

  if (xi_given) {
    fq0 <- box_cox(-log(q0), lambda = -xi)
    fq2 <- box_cox(-log(q2), lambda = -xi)
    sigma <- (x[3] - x[1]) / (fq0 - fq2)
    mu <- x[1] + sigma * fq0
  } else {
    r_hat <- (x[3] - x[2]) / (x[2] - x[1])
    xi <- -log(r_hat) / log(a)
    if (sigma_given) {
      if (xi < 0) {
        mu <- x[3] + sigma * box_cox(-log(q2), lambda = -xi)
      } else {
        mu <- x[1] + sigma * box_cox(-log(q0), lambda = -xi)
      }
    } else if (mu_given) {
      if (xi < 0) {
        sigma <- (mu - x[3]) / box_cox(-log(q2), lambda = -xi)
      } else {
        sigma <- (x[1] - mu) / box_cox(-log(q0), lambda = xi)
      }
    }
  }
  val <- c(mu, sigma, xi)
  names(val) <- c("mu", "sigma", "xi")
  return(val)
}
