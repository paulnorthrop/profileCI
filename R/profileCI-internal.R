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
profile_ci <- function(object, negated_loglik_fn, which = 1, level, mle, inc,
                       epsilon, ...) {

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
    opt <- try(stats::optim(sol, profiling_fn, par_which = par_which, ...),
               silent = TRUE)
    # If optim errors then set different initial values and try again
    if (inherits(opt, "try-error")) {
      init <- initial_mvn(object = object, x = which, x_value = par_which)
      opt <- try(stats::optim(init, profiling_fn, par_which = par_which, ...),
                 silent = TRUE)
    }
    if (inherits(opt, "try-error")) {
      return(list(optim_error = attr(opt, "condition")))
    }
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
    opt <- try(stats::optim(sol, profiling_fn, par_which = par_which, ...),
               silent = TRUE)
    # If optim errors then set different initial values and try again
    if (inherits(opt, "try-error")) {
      init <- initial_mvn(object = object, x = which, x_value = par_which)
      opt <- try(stats::optim(init, profiling_fn, par_which = par_which, ...),
                 silent = TRUE)
    }
    if (inherits(opt, "try-error")) {
      return(list(optim_error = attr(opt, "condition")))
    }
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

  # If epsilon = 0 then use linear interpolation
  # If epsilon < 0 then use quadratic interpolation
  # If epsilon > 0 then use quadratic interpolation and then itp::itp()

  up_lim <- x1up + (conf_line - y1up) * (x2up - x1up) / (y2up - y1up)
  low_lim <- x1low + (conf_line - y1low) * (x2low - x1low) / (y2low - y1low)
  if (epsilon != 0) {
    # Calculate the values of the profile log-likelihood at these limits and
    # use the 3 points (the bracketing points and this new point) to estimate
    # the confidence limits by quadratic interpolation

    # Upper
    opt <- try(stats::optim(sol_upper, profiling_fn,
                            par_which = up_lim, ...), silent = TRUE)
    # If optim errors then set different initial values and try again
    if (inherits(opt, "try-error")) {
      init <- initial_mvn(object = object, x = which, x_value = up_lim)
      opt <- try(stats::optim(init, profiling_fn, par_which = up_lim, ...),
                 silent = TRUE)
    }
    up_new <- -opt$value
    temp <- lagrangianInterpolation(c(y1up, up_new, y2up),
                                    c(x1up, up_lim, x2up))
    save_up_lim <- up_lim
    up_lim <- temp(conf_line)

    # Lower
    opt <- try(stats::optim(sol_lower, profiling_fn,
                            par_which = low_lim, ...), silent = TRUE)
    # If optim errors then set different initial values and try again
    if (inherits(opt, "try-error")) {
      init <- initial_mvn(object = object, x = which, x_value = low_lim)
      opt <- try(stats::optim(init, profiling_fn, par_which = low_lim, ...),
                 silent = TRUE)
    }
    low_new <- -opt$value
    temp <- lagrangianInterpolation(c(y1low, low_new, y2low),
                                    c(x1low, low_lim, x2low))
    save_low_lim <- low_lim
    low_lim <- temp(conf_line)

    # If epsilon > 0 then use itp::itp(), creating a new bracket from 2 of the
    # 3 points available and set an initial estimate
    if (epsilon > 0) {
      itp_function <- function(x, par, ...) {
        opt <- try(stats::optim(par = par, fn = profiling_fn, par_which = x,
                                ...), silent = TRUE)
        if (inherits(opt, "try-error")) {
          init <- initial_mvn(object = object, x = which, x_value = x)
          opt <- try(stats::optim(par = init, fn = profiling_fn, par_which = x,
                                  ...), silent = TRUE)
        }
        val <- -opt$value - conf_line
        attr(val, "parameters") <- opt$par
        return(val)
      }
      # Find the upper limit of the confidence interval
      if (up_new > 0) {
        interval <- c(x1up, save_up_lim)
      } else {
        interval <- c(save_up_lim, x2up)
      }
      upper <- itp::itp(f = itp_function, interval = interval,
                        par = sol_upper,
                        f.a = y1up - conf_line, f.b = y2up - conf_line,
                        epsilon = epsilon, ...)
      up_lim <- upper$root
      # Find the lower limit of the confidence interval
      if (low_new > 0) {
        interval <- c(x1low, save_low_lim)
      } else {
        interval <- c(save_low_lim, x2low)
      }
      lower <- itp::itp(f = itp_function, interval = interval,
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
                    prof_lik[(loc_lower + 1):loc_upper],
                    upper$f.root + conf_line, prof_lik[(loc_upper + 1):n])
      # Save the parameter values that apply to the solutions from itp::itp()
      lower_pars <- numeric(n_pars)
      lower_pars[which] <- lower$root
      lower_pars[-which] <- attr(lower$f.root, "parameters")
      names(lower_pars) <- names(mle)
      upper_pars <- numeric(n_pars)
      upper_pars[which] <- upper$root
      upper_pars[-which] <- attr(upper$f.root, "parameters")
      names(upper_pars) <- names(mle)
    }
  }

  par_prof <- c(lower = low_lim, mle_which, upper = up_lim)
  return(list(par_prof = par_prof, crit = conf_line,
              for_plot = cbind(par_values = par_values,
                               prof_loglik = prof_lik),
              lower_pars = lower_pars, upper_pars = upper_pars))
}

# ========================== Faster profiling function ====================== #

#' @keywords internal
#' @rdname profileCI-internal
faster_profile_ci <- function(object, negated_loglik_fn, which = 1, which_name,
                              level, mle, ci_sym_mat, inc, epsilon, ...) {

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
  # initial_mvn() to set initial estimates of the parameters other than the
  # parameter which, based on the estimated approximate multivariate normal
  # distribution of the MLE

  par_which <- ci_sym_mat[which_name, 2]
  init <- initial_mvn(object = object, x = which, x_value = par_which)

  # Calculate the profile log-likelihood at the initial values
  # If this is greater than conf_line then we search upwards
  # If this is smaller than conf_line then we search downwards

  ii <- 2
  opt <- try(stats::optim(init, profiling_fn, par_which = par_which, ...),
             silent = TRUE)
  # If optim errors then set different initial values and try again
  if (inherits(opt, "try-error")) {
    init <- mle[-which]
    opt <- try(stats::optim(init, profiling_fn, par_which = par_which, ...),
               silent = TRUE)
  }
  if (inherits(opt, "try-error")) {
    return(list(optim_error = attr(opt, "condition")))
  }
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
    opt <- try(stats::optim(sol, profiling_fn, par_which = par_which, ...),
               silent = TRUE)
    # If optim errors then reset the initial values and try again
    if (inherits(opt, "try-error")) {
      init <- initial_mvn(object = object, x = which, x_value = par_which)
      opt <- try(stats::optim(init, profiling_fn, par_which = par_which, ...),
                 silent = TRUE)
    }
    if (inherits(opt, "try-error")) {
      return(list(optim_error = attr(opt, "condition")))
    }
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

  # We start from the lower limit of the symmetric confidence interval, using
  # initial_mvn() to set initial estimates of the parameters other than the
  # parameter which, based on the estimated approximate multivariate normal
  # distribution of the MLE

  par_which <- ci_sym_mat[which_name, 1]
  init <- initial_mvn(object = object, x = which, x_value = par_which)

  # Calculate the profile log-likelihood at the initial values
  # If this is greater than conf_line then we search downwards
  # If this is smaller than conf_line then we search upwards

  ii <- 2
  opt <- try(stats::optim(init, profiling_fn, par_which = par_which, ...),
             silent = TRUE)
  # If optim errors then set different initial values and try again
  if (inherits(opt, "try-error")) {
    init <- mle[-which]
    opt <- try(stats::optim(init, profiling_fn, par_which = par_which, ...),
               silent = TRUE)
  }
  if (inherits(opt, "try-error")) {
    return(list(optim_error = attr(opt, "condition")))
  }
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
    opt <- try(stats::optim(sol, profiling_fn, par_which = par_which, ...),
               silent = TRUE)
    # If optim errors then reset the initial values and try again
    if (inherits(opt, "try-error")) {
      init <- initial_mvn(object = object, x = which, x_value = par_which)
      opt <- try(stats::optim(init, profiling_fn, par_which = par_which, ...),
                 silent = TRUE)
    }
    if (inherits(opt, "try-error")) {
      return(list(optim_error = attr(opt, "condition")))
    }
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

  # If epsilon = 0 then use linear interpolation
  # If epsilon < 0 then use quadratic interpolation
  # If epsilon > 0 then use quadratic interpolation and then itp::itp()

  lower_pars <- NULL
  upper_pars <- NULL
  up_lim <- x1up + (conf_line - y1up) * (x2up - x1up) / (y2up - y1up)
  low_lim <- x1low + (conf_line - y1low) * (x2low - x1low) / (y2low - y1low)
  if (epsilon != 0) {
    # Calculate the values of the profile log-likelihood at these limits and
    # use the 3 points (the bracketing points and this new point) to estimate
    # the confidence limits by quadratic interpolation

    # Upper
    opt <- try(stats::optim(sol_upper, profiling_fn,
                            par_which = up_lim, ...), silent = TRUE)
    # If optim errors then reset the initial values and try again
    if (inherits(opt, "try-error")) {
      init <- initial_mvn(object = object, x = which, x_value = up_lim)
      opt <- try(stats::optim(init, profiling_fn, par_which = up_lim, ...),
                 silent = TRUE)
    }
    up_new <- -opt$value
    temp <- lagrangianInterpolation(c(y1up, up_new, y2up),
                                    c(x1up, up_lim, x2up))
    save_up_lim <- up_lim
    up_lim <- temp(conf_line)

    # Lower
    opt <- try(stats::optim(sol_lower, profiling_fn,
                            par_which = low_lim, ...), silent = TRUE)
    # If optim errors then reset the initial values and try again
    if (inherits(opt, "try-error")) {
      init <- initial_mvn(object = object, x = which, x_value = low_lim)
      opt <- try(stats::optim(init, profiling_fn, par_which = low_lim, ...),
                 silent = TRUE)
    }
    low_new <- -opt$value
    temp <- lagrangianInterpolation(c(y1low, low_new, y2low),
                                    c(x1low, low_lim, x2low))
    save_low_lim <- low_lim
    low_lim <- temp(conf_line)

    # Add the approximate solutions to par_values and prof_lik
    # Note that only the extreme ends of these vectors have prof_lik values
    # below conf_line
    # Only do this if epsilon < 0 to avoid messing up the ordering of points
    # when we do something similar below for the epsilon > 0 case
    if (epsilon < 0) {
      n <- length(par_values)
      par_values <- c(par_values[1:loc_lower], low_lim,
                      par_values[(loc_lower + 1):loc_upper], up_lim,
                      par_values[(loc_upper + 1):n])
      prof_lik <- c(prof_lik[1:loc_lower], low_new,
                    prof_lik[(loc_lower + 1):loc_upper],
                    up_new, prof_lik[(loc_upper + 1):n])
    }

    # If epsilon > 0 then use itp::itp(), creating a new bracket from 2 of the
    # 3 points available and set an initial estimate
    if (epsilon > 0) {
      itp_function <- function(x, par, ...) {
        opt <- try(stats::optim(par = par, fn = profiling_fn, par_which = x,
                                ...), silent = TRUE)
        if (inherits(opt, "try-error")) {
          init <- initial_mvn(object = object, x = which, x_value = x)
          opt <- try(stats::optim(par = init, fn = profiling_fn, par_which = x,
                                  ...), silent = TRUE)
        }
        val <- -opt$value - conf_line
        attr(val, "parameters") <- opt$par
        return(val)
      }
      # Find the upper limit of the confidence interval
      if (up_new > 0) {
        interval <- c(x1up, save_up_lim)
      } else {
        interval <- c(save_up_lim, x2up)
      }
      upper <- itp::itp(f = itp_function, interval = interval,
                        par = sol_upper,
                        f.a = y1up - conf_line, f.b = y2up - conf_line,
                        epsilon = epsilon, ...)
      up_lim <- upper$root
      # Find the lower limit of the confidence interval
      if (low_new > 0) {
        interval <- c(x1low, save_low_lim)
      } else {
        interval <- c(save_low_lim, x2low)
      }
      lower <- itp::itp(f = itp_function, interval = interval,
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
                    prof_lik[(loc_lower + 1):loc_upper],
                    upper$f.root + conf_line, prof_lik[(loc_upper + 1):n])
      # Save the parameter values that apply to the solutions from itp::itp()
      lower_pars <- numeric(n_pars)
      lower_pars[which] <- lower$root
      lower_pars[-which] <- attr(lower$f.root, "parameters")
      names(lower_pars) <- names(mle)
      upper_pars <- numeric(n_pars)
      upper_pars[which] <- upper$root
      upper_pars[-which] <- attr(upper$f.root, "parameters")
      names(upper_pars) <- names(mle)
    }
  }

  par_prof <- c(lower = low_lim, mle_which, upper = up_lim)
  return(list(par_prof = par_prof, crit = conf_line,
              for_plot = cbind(par_values = par_values,
                               prof_loglik = prof_lik),
              lower_pars = lower_pars, upper_pars = upper_pars))
}

#' @keywords internal
#' @rdname profileCI-internal
lagrangianInterpolation <- function(x0, y0) {
  f <- function(x) {
    sum(y0 * sapply(seq_along(x0), \(j) {
      prod(x - x0[-j])/prod(x0[j] - x0[-j])
    }))
  }
  return(Vectorize(f, "x"))
}

#' @keywords internal
#' @rdname profileCI-internal
initial_mvn <- function(object, x, x_value) {
  # MLE and variance-covariance matrix
  mle <- coef(object)
  cov <- vcov(object)
  # Calculate the conditional mean of y[-x] | x = x_value
  y <- (1:length(mle))[-x]
  sigma12 <- cov[y, x, drop = FALSE]
  sigma22 <- cov[x, x]
  cond_mean <- c(mle[y] + sigma12 %*% solve(sigma22) %*% (x_value - mle[x]))
  return(cond_mean)
}
