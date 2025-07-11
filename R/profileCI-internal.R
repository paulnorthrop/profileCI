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
                       epsilon, optim_args, mult, flat, lb, ub, ...) {
  # To avoid potential issues with passing arguments to the negated
  # log-likelihood function via ... we save these arguments in loglik_args now
  # and use them inside profiling_fn()
  loglik_args <- list(...)

  # The number of parameters
  n_pars <- length(mle)

  # The -log-likelihood to profile over parameters in par other than par[which]
  profiling_fn <- function(par, par_which) {
    # Place par_which in position which and the other parameters around it
    parameters <- rep_len(NA, length(par) + 1)
    parameters[which] <- par_which
    parameters[-which] <- par
    # Return the negated log-likelihood
    negated_loglik_fn_args <- c(list(parm = parameters), loglik_args)
    val <- do.call(negated_loglik_fn, negated_loglik_fn_args)
    return(val)
  }

  # The maximised log-likelihood and the MLE for the parameter of interest
  mle_args <- c(list(parm = mle), loglik_args)
  max_loglik <- -do.call(negated_loglik_fn, mle_args)
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

  # Infer 1 SE from inc and mult for later use
  one_se <- 100 * inc / mult

  ### Upper tail ...

  par_which <- mle_which
  my_val <- max_loglik
  ii <- 1
  sol <- mle[-which]
  # Add a check for flatness of the profile log-likelihood
  flat_upper <- FALSE
  # Add a check for hitting the upper bound ub
  hit_ub <- FALSE
  while (my_val > conf_line && !flat_upper && !hit_ub){
    par_which <- par_which + inc
    o_args <- list(par = sol, fn = profiling_fn, par_which = par_which)
    opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
    # If optim errors then set different initial values and try again
    if (inherits(opt, "try-error")) {
      init <- initial_mvn(object = object, x = which, x_value = par_which)
      o_args <- list(par = init, fn = profiling_fn, par_which = par_which)
      opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
    }
    if (inherits(opt, "try-error")) {
      return(list(optim_error = attr(opt, "condition")))
    }
    sol <- opt$par
    ii <- ii + 1
    x2[ii] <- par_which
    v2[ii] <- -opt$value
    # Check for flatness
    # Use one_se_away to avoid stopping near the MLE
    one_se_away <- abs(par_which - mle_which) > one_se
    delta_loglik <- 100 * (my_val + opt$value) / mult
    if (one_se_away && delta_loglik > 0 && delta_loglik < flat) {
      flat_upper <- TRUE
    }
    if (par_which >= ub) {
      hit_ub <- TRUE
    }
    # Save the current value of the profile log-likelihood
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
  # Add a check for flatness of the profile log-likelihood
  flat_lower <- FALSE
  # Add a check for hitting the lower bound lb
  hit_lb <- FALSE
  while (my_val > conf_line && !flat_lower && !hit_lb){
    par_which <- par_which - inc
    o_args <- list(par = sol, fn = profiling_fn, par_which = par_which)
    opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
    # If optim errors then set different initial values and try again
    if (inherits(opt, "try-error")) {
      init <- initial_mvn(object = object, x = which, x_value = par_which)
      o_args <- list(par = init, fn = profiling_fn, par_which = par_which)
      opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
    }
    if (inherits(opt, "try-error")) {
      return(list(optim_error = attr(opt, "condition")))
    }
    sol <- opt$par
    ii <- ii + 1
    x1[ii] <- par_which
    v1[ii] <- -opt$value
    # Check for flatness
    # Use one_se_away to avoid stopping near the MLE
    one_se_away <- abs(par_which - mle_which) > one_se
    delta_loglik <- 100 * (my_val + opt$value) / mult
    if (one_se_away && delta_loglik > 0 && delta_loglik < flat) {
      flat_lower <- TRUE
    }
    if (par_which <= lb) {
      hit_lb <- TRUE
    }
    # Save the current value of the profile log-likelihood
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
  if (flat_upper) {
    up_lim <- Inf
  } else if (hit_ub) {
    up_lim <- NA
  } else {
    loc_upper <- which(temp == -1)
    x1up <- par_values[loc_upper]
    x2up <- par_values[loc_upper + 1]
    y1up <- prof_lik[loc_upper]
    y2up <- prof_lik[loc_upper + 1]
    up_lim <- x1up + (conf_line - y1up) * (x2up - x1up) / (y2up - y1up)
  }
  # Find the lower limit of the confidence interval
  if (flat_lower) {
    low_lim <- -Inf
  } else if (hit_lb) {
    low_lim <- NA
  } else {
    loc_lower <- which(temp == 1)
    x1low <- par_values[loc_lower]
    x2low <- par_values[loc_lower + 1]
    y1low <- prof_lik[loc_lower]
    y2low <- prof_lik[loc_lower + 1]
    low_lim <- x1low + (conf_line - y1low) * (x2low - x1low) / (y2low - y1low)
  }

  # If epsilon = 0 then use linear interpolation
  # If epsilon < 0 then use quadratic interpolation
  # If epsilon > 0 then use quadratic interpolation and then itp::itp()

  if (epsilon != 0) {
    # Calculate the values of the profile log-likelihood at these limits and
    # use the 3 points (the bracketing points and this new point) to estimate
    # the confidence limits by quadratic interpolation

    # Upper
    if (flat_upper) {
      up_lim <- Inf
    } else if (hit_ub) {
      up_lim <- NA
    } else {
      o_args <- list(par = sol_upper, fn = profiling_fn, par_which = up_lim)
      opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
      # If optim errors then set different initial values and try again
      if (inherits(opt, "try-error")) {
        init <- initial_mvn(object = object, x = which, x_value = up_lim)
        o_args <- list(par = init, fn = profiling_fn, par_which = up_lim)
        opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
      }
      up_new <- -opt$value
      temp <- lagrangianInterpolation(c(y1up, up_new, y2up),
                                      c(x1up, up_lim, x2up))
      save_up_lim <- up_lim
      up_lim <- temp(conf_line)
    }

    # Lower
    if (flat_lower) {
      low_lim <- -Inf
    } else if (hit_lb) {
      low_lim <- NA
    } else {
      o_args <- list(par = sol_lower, fn = profiling_fn, par_which = low_lim)
      opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
      # If optim errors then set different initial values and try again
      if (inherits(opt, "try-error")) {
        init <- initial_mvn(object = object, x = which, x_value = low_lim)
        o_args <- list(par = init, fn = profiling_fn, par_which = low_lim)
        opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
      }
      low_new <- -opt$value
      temp <- lagrangianInterpolation(c(y1low, low_new, y2low),
                                      c(x1low, low_lim, x2low))
      save_low_lim <- low_lim
      low_lim <- temp(conf_line)
    }

    # If epsilon > 0 then use itp::itp(), creating a new bracket from 2 of the
    # 3 points available and set an initial estimate
    if (epsilon > 0) {
      itp_function <- function(val_par_which, par) {
        o_args <- list(par = par, fn = profiling_fn, par_which = val_par_which)
        opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
        if (inherits(opt, "try-error")) {
          init <- initial_mvn(object = object, x = which,
                              x_value = val_par_which)
          o_args <- list(par = init, fn = profiling_fn,
                         par_which = val_par_which)
          opt <- try(do.call(stats::optim, c(o_args, optim_args)),
                     silent = TRUE)
        }
        val <- -opt$value - conf_line
        attr(val, "parameters") <- opt$par
        return(val)
      }
      if (flat_upper) {
        up_lim <- Inf
      } else if (hit_ub) {
        up_lim <- NA
      } else {
        # Find the upper limit of the confidence interval
        if (up_new > 0) {
          interval <- c(x1up, save_up_lim)
        } else {
          interval <- c(save_up_lim, x2up)
        }
        upper <- itp::itp(f = itp_function, interval = interval,
                          par = sol_upper,
                          f.a = y1up - conf_line, f.b = y2up - conf_line,
                          epsilon = epsilon)
        up_lim <- upper$root
      }
      if (flat_lower) {
        low_lim <- -Inf
      } else if (hit_lb) {
        low_lim <- NA
      } else {
        # Find the lower limit of the confidence interval
        if (low_new > 0) {
          interval <- c(x1low, save_low_lim)
        } else {
          interval <- c(save_low_lim, x2low)
        }
        lower <- itp::itp(f = itp_function, interval = interval,
                          par = sol_lower,
                          f.a = y1low - conf_line, f.b = y2low - conf_line,
                          epsilon = epsilon)
        low_lim <- lower$root
      }

      # Add the approximate solutions to par_values and prof_lik
      # Note that only the extreme ends of these vectors have prof_lik values
      # below conf_line

      n <- length(par_values)
      if (flat_upper || hit_ub) {
        add_upper_lim <- NULL
        add_upper_prof <- NULL
        loc_upper <- n - 1
      } else {
        add_upper_lim <- up_lim
        add_upper_prof <- upper$f.root + conf_line
      }
      if (flat_lower || hit_lb) {
        add_lower_lim <- NULL
        add_lower_prof <- NULL
        loc_lower <- 1
      } else {
        add_lower_lim <- low_lim
        add_lower_prof <- lower$f.root + conf_line
      }
      par_values <- c(par_values[1:loc_lower], add_lower_lim,
                      par_values[(loc_lower + 1):loc_upper], add_upper_lim,
                      par_values[(loc_upper + 1):n])
      prof_lik <- c(prof_lik[1:loc_lower], add_lower_prof,
                    prof_lik[(loc_lower + 1):loc_upper],
                    add_upper_prof, prof_lik[(loc_upper + 1):n])
      # Save the parameter values that apply to the solutions from itp::itp()
      if (!flat_lower && !hit_lb) {
        lower_pars <- numeric(n_pars)
        lower_pars[which] <- lower$root
        lower_pars[-which] <- attr(lower$f.root, "parameters")
        names(lower_pars) <- names(mle)
      }
      if (!flat_upper && !hit_ub) {
        upper_pars <- numeric(n_pars)
        upper_pars[which] <- upper$root
        upper_pars[-which] <- attr(upper$f.root, "parameters")
        names(upper_pars) <- names(mle)
      }
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
                              level, mle, ci_sym_mat, inc, epsilon, optim_args,
                              mult, flat, lb, ub, ...) {
  # To avoid potential issues with passing arguments to the negated
  # log-likelihood function via ... we save these arguments in loglik_args now
  # and use them inside profiling_fn()
  loglik_args <- list(...)

  # The number of parameters
  n_pars <- length(mle)

  # The -log-likelihood to profile over parameters in par other than par[which]
  profiling_fn <- function(par, par_which) {
    # Place par_which in position which and the other parameters around it
    parameters <- rep_len(NA, length(par) + 1)
    parameters[which] <- par_which
    parameters[-which] <- par
    # Return the negated log-likelihood
    negated_loglik_fn_args <- c(list(parm = parameters), loglik_args)
    val <- do.call(negated_loglik_fn, negated_loglik_fn_args)
    return(val)
  }

  # The maximised log-likelihood and the MLE for the parameter of interest
  mle_args <- c(list(parm = mle), loglik_args)
  max_loglik <- -do.call(negated_loglik_fn, mle_args)
  mle_which <- mle[which]
  # The horizontal line that determines the confidence limits
  conf_line <- max_loglik - 0.5 * stats::qchisq(level, 1)
  # Vectors to store values of the parameters (x1 and x2) and the values of
  # the profile log-likelihood (v1 and v2)
  v1 <- v2 <- x1 <- x2 <- NULL
  x2[1] <- x1[1] <- mle[which]
  v2[1] <- v1[1] <- max_loglik

  # Infer 1 SE from inc and mult for later use
  one_se <- 100 * inc / mult

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
  o_args <- list(par = init, fn = profiling_fn, par_which = par_which)
  opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
  # If optim errors then set different initial values and try again
  if (inherits(opt, "try-error")) {
    init <- mle[-which]
    o_args <- list(par = init, fn = profiling_fn, par_which = par_which)
    opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
  }
  if (inherits(opt, "try-error")) {
    return(list(optim_error = attr(opt, "condition")))
  }
  my_val <- -opt$value
  x2[ii] <- par_which
  v2[ii] <- my_val

  # If my_val < conf_line then the profile log-likelihood has dropped below
  # the required level in one (big) step. We use linear interpolation between
  # this point and the MLE to move back (hopefully, just) above the level.
  # This should work because the profile log-likelihood should be convex.
  if (my_val < conf_line) {
    searched_upwards <- FALSE
    while_condition <- function(my_val) {
      return(my_val < conf_line)
    }
    temp <- lagrangianInterpolation(c(v2[1], v2[2]), c(x2[1], x2[2]))
    delta <- -(x2[2] - temp(conf_line))
  } else {
    searched_upwards <- TRUE
    while_condition <- function(my_val) {
      return(my_val > conf_line)
    }
    delta <- inc
  }
  sol <- opt$par

  # Add a check for flatness of the profile log-likelihood
  flat_upper <- FALSE
  # Add a check for hitting the upper bound ub
  hit_ub <- FALSE
  while (while_condition(my_val) && !flat_upper && !hit_ub){
    par_which <- par_which + delta
    o_args <- list(par = sol, fn = profiling_fn, par_which = par_which)
    opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
    # If optim errors then reset the initial values and try again
    if (inherits(opt, "try-error")) {
      init <- initial_mvn(object = object, x = which, x_value = par_which)
      o_args <- list(par = init, fn = profiling_fn, par_which = par_which)
      opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
    }
    if (inherits(opt, "try-error")) {
      return(list(optim_error = attr(opt, "condition")))
    }
    sol <- opt$par
    ii <- ii + 1
    x2[ii] <- par_which
    v2[ii] <- -opt$value
    # Check for flatness
    # Use three_se_away to avoid stopping near the MLE
    three_se_away <- abs(par_which - mle_which) > 3 * one_se
    delta_loglik <- 100 * (my_val + opt$value) / mult
    if (three_se_away && delta_loglik > 0 && delta_loglik < flat) {
      flat_upper <- TRUE
    }
    if (par_which >= ub) {
      hit_ub <- TRUE
    }
    # Save the current value of the profile log-likelihood
    my_val <- v2[ii]
  }
  # If we searched downwards then reorder the results
  if (!searched_upwards) {
    x2[-1] <- rev(x2[-1])
    v2[-1] <- rev(v2[-1])
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
  o_args <- list(par = init, fn = profiling_fn, par_which = par_which)
  opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
  # If optim errors then set different initial values and try again
  if (inherits(opt, "try-error")) {
    init <- mle[-which]
    o_args <- list(par = init, fn = profiling_fn, par_which = par_which)
    opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
  }
  if (inherits(opt, "try-error")) {
    return(list(optim_error = attr(opt, "condition")))
  }
  my_val <- -opt$value
  x1[ii] <- par_which
  v1[ii] <- my_val

  # If my_val < conf_line then the profile log-likelihood has dropped below
  # the required level in one (big) step. We use linear interpolation between
  # this point and the MLE to move back (hopefully, just) above the level.
  # This should work because the profile log-likelihood should be convex.
  if (my_val < conf_line) {
    searched_downwards <- FALSE
    while_condition <- function(my_val) {
      return(my_val < conf_line)
    }
    temp <- lagrangianInterpolation(c(v1[1], v1[2]), c(x1[1], x1[2]))
    delta <- -(temp(conf_line) - x1[2])
  } else {
    searched_downwards <- TRUE
    while_condition <- function(my_val) {
      return(my_val > conf_line)
    }
    delta <- inc
  }
  sol <- opt$par

  # Add a check for flatness of the profile log-likelihood
  flat_lower <- FALSE
  # Add a check for hitting the lower bound lb
  hit_lb <- FALSE
  while (while_condition(my_val) && !flat_lower && !hit_lb){
    par_which <- par_which - delta
    o_args <- list(par = sol, fn = profiling_fn, par_which = par_which)
    opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
    # If optim errors then reset the initial values and try again
    if (inherits(opt, "try-error")) {
      init <- initial_mvn(object = object, x = which, x_value = par_which)
      o_args <- list(par = init, fn = profiling_fn, par_which = par_which)
      opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
    }
    if (inherits(opt, "try-error")) {
      return(list(optim_error = attr(opt, "condition")))
    }
    sol <- opt$par
    ii <- ii + 1
    x1[ii] <- par_which
    v1[ii] <- -opt$value
    # Check for flatness
    # Use three_se_away to avoid stopping near the MLE
    three_se_away <- abs(par_which - mle_which) > 3 * one_se
    delta_loglik <- 100 * (my_val + opt$value) / mult
    if (three_se_away && delta_loglik > 0 && delta_loglik < flat) {
      flat_lower <- TRUE
    }
    if (par_which <= lb) {
      hit_lb <- TRUE
    }
    # Save the current value of the profile log-likelihood
    my_val <- v1[ii]
  }
  # If we searched upwards then reorder the results
  if (!searched_downwards) {
    x1[-1] <- rev(x1[-1])
    v1[-1] <- rev(v1[-1])
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
  if (flat_upper) {
    up_lim <- Inf
  } else if (hit_ub) {
    up_lim <- NA
  } else {
    loc_upper <- which(temp == -1)
    x1up <- par_values[loc_upper]
    x2up <- par_values[loc_upper + 1]
    y1up <- prof_lik[loc_upper]
    y2up <- prof_lik[loc_upper + 1]
    up_lim <- x1up + (conf_line - y1up) * (x2up - x1up) / (y2up - y1up)

  }
  # Find the lower limit of the confidence interval
  if (flat_lower) {
    low_lim <- -Inf
  } else if (hit_lb) {
    low_lim <- NA
  } else {
    loc_lower <- which(temp == 1)
    x1low <- par_values[loc_lower]
    x2low <- par_values[loc_lower + 1]
    y1low <- prof_lik[loc_lower]
    y2low <- prof_lik[loc_lower + 1]
    low_lim <- x1low + (conf_line - y1low) * (x2low - x1low) / (y2low - y1low)
  }

  # If epsilon = 0 then use linear interpolation
  # If epsilon < 0 then use quadratic interpolation
  # If epsilon > 0 then use quadratic interpolation and then itp::itp()

  if (epsilon != 0) {
    # Calculate the values of the profile log-likelihood at these limits and
    # use the 3 points (the bracketing points and this new point) to estimate
    # the confidence limits by quadratic interpolation

    # Upper
    if (flat_upper) {
      up_lim <- Inf
    } else if (hit_ub) {
      up_lim <- NA
    } else {
      o_args <- list(par = sol_upper, fn = profiling_fn, par_which = up_lim)
      opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
      # If optim errors then reset the initial values and try again
      if (inherits(opt, "try-error")) {
        init <- initial_mvn(object = object, x = which, x_value = up_lim)
        o_args <- list(par = init, fn = profiling_fn, par_which = up_lim)
        opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
      }
      up_new <- -opt$value
      temp <- lagrangianInterpolation(c(y1up, up_new, y2up),
                                      c(x1up, up_lim, x2up))
      save_up_lim <- up_lim
      up_lim <- temp(conf_line)
    }

    # Lower
    if (flat_lower) {
      low_lim <- -Inf
    } else if (hit_lb) {
      low_lim <- NA
    } else {
      o_args <- list(par = sol_lower, fn = profiling_fn, par_which = low_lim)
      opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
      # If optim errors then reset the initial values and try again
      if (inherits(opt, "try-error")) {
        init <- initial_mvn(object = object, x = which, x_value = low_lim)
        o_args <- list(par = init, fn = profiling_fn, par_which = low_lim)
        opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
      }
      low_new <- -opt$value
      temp <- lagrangianInterpolation(c(y1low, low_new, y2low),
                                      c(x1low, low_lim, x2low))
      save_low_lim <- low_lim
      low_lim <- temp(conf_line)
    }

    # Add the approximate solutions to par_values and prof_lik
    # Note that only the extreme ends of these vectors have prof_lik values
    # below conf_line
    # Only do this if epsilon < 0 to avoid messing up the ordering of points
    # when we do something similar below for the epsilon > 0 case
    if (epsilon < 0) {
      n <- length(par_values)
      if (flat_upper || hit_ub) {
        add_upper_lim <- NULL
        add_upper_prof <- NULL
        loc_upper <- n - 1
      } else {
        add_upper_lim <- up_lim
        add_upper_prof <- up_new
      }
      if (flat_lower || hit_lb) {
        add_lower_lim <- NULL
        add_lower_prof <- NULL
        loc_lower <- 1
      } else {
        add_lower_lim <- low_lim
        add_lower_prof <- low_new
      }
      par_values <- c(par_values[1:loc_lower], add_lower_lim,
                      par_values[(loc_lower + 1):loc_upper], add_upper_lim,
                      par_values[(loc_upper + 1):n])
      prof_lik <- c(prof_lik[1:loc_lower], add_lower_prof,
                    prof_lik[(loc_lower + 1):loc_upper],
                    add_upper_prof, prof_lik[(loc_upper + 1):n])
    }

    # If epsilon > 0 then use itp::itp(), creating a new bracket from 2 of the
    # 3 points available and set an initial estimate
    if (epsilon > 0) {
      itp_function <- function(val_par_which, par) {
        o_args <- list(par = par, fn = profiling_fn, par_which = val_par_which)
        opt <- try(do.call(stats::optim, c(o_args, optim_args)), silent = TRUE)
        if (inherits(opt, "try-error")) {
          init <- initial_mvn(object = object, x = which,
                              x_value = val_par_which)
          o_args <- list(par = init, fn = profiling_fn,
                         par_which = val_par_which)
          opt <- try(do.call(stats::optim, c(o_args, optim_args)),
                     silent = TRUE)
        }
        val <- -opt$value - conf_line
        attr(val, "parameters") <- opt$par
        return(val)
      }
      if (flat_upper) {
        up_lim <- Inf
      } else if (hit_ub) {
        up_lim <- NA
      } else {
        # Find the upper limit of the confidence interval
        if (up_new > 0) {
          interval <- c(x1up, save_up_lim)
        } else {
          interval <- c(save_up_lim, x2up)
        }
        upper <- itp::itp(f = itp_function, interval = interval,
                          par = sol_upper,
                          f.a = y1up - conf_line, f.b = y2up - conf_line,
                          epsilon = epsilon)
        up_lim <- upper$root
      }
      if (flat_lower) {
        low_lim <- -Inf
      } else if (hit_lb) {
        low_lim <- NA
      } else {
        # Find the lower limit of the confidence interval
        if (low_new > 0) {
          interval <- c(x1low, save_low_lim)
        } else {
          interval <- c(save_low_lim, x2low)
        }
        lower <- itp::itp(f = itp_function, interval = interval,
                          par = sol_lower,
                          f.a = y1low - conf_line, f.b = y2low - conf_line,
                          epsilon = epsilon)
        low_lim <- lower$root
      }
      # Add the approximate solutions to par_values and prof_lik
      # Note that only the extreme ends of these vectors have prof_lik values
      # below conf_line

      n <- length(par_values)
      if (flat_upper || hit_ub) {
        add_upper_lim <- NULL
        add_upper_prof <- NULL
        loc_upper <- n - 1
      } else {
        add_upper_lim <- up_lim
        add_upper_prof <- upper$f.root + conf_line
      }
      if (flat_lower || hit_lb) {
        add_lower_lim <- NULL
        add_lower_prof <- NULL
        loc_lower <- 1
      } else {
        add_lower_lim <- low_lim
        add_lower_prof <- lower$f.root + conf_line
      }
      par_values <- c(par_values[1:loc_lower], add_lower_lim,
                      par_values[(loc_lower + 1):loc_upper], add_upper_lim,
                      par_values[(loc_upper + 1):n])
      prof_lik <- c(prof_lik[1:loc_lower], add_lower_prof,
                    prof_lik[(loc_lower + 1):loc_upper],
                    add_upper_prof, prof_lik[(loc_upper + 1):n])

      # Save the parameter values that apply to the solutions from itp::itp()
      if (!flat_lower && !hit_lb) {
        lower_pars <- numeric(n_pars)
        lower_pars[which] <- lower$root
        lower_pars[-which] <- attr(lower$f.root, "parameters")
        names(lower_pars) <- names(mle)
      }
      if (!flat_upper && !hit_ub) {
        upper_pars <- numeric(n_pars)
        upper_pars[which] <- upper$root
        upper_pars[-which] <- attr(upper$f.root, "parameters")
        names(upper_pars) <- names(mle)
      }
    }
  }

  par_prof <- c(lower = low_lim, mle_which, upper = up_lim)
  return(list(par_prof = par_prof, crit = conf_line,
              for_plot = cbind(par_values = par_values,
                               prof_loglik = prof_lik),
              lower_pars = lower_pars, upper_pars = upper_pars))
}

# ========================= Lagrangian interpolation ======================== #

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

# ======================== Initial estimates from vcov ====================== #

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
