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
profile_ci <- function(negated_loglik_fn, which = 1, level, mle, inc, ...) {

  # The -log-likelihood to profile over parameters in par other than par[which]
  profiling_fn <- function(par, ...) {
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
    opt <- stats::optim(sol, profiling_fn, ...)
    sol <- opt$par
    ii <- ii + 1
    x2[ii] <- par_which
    v2[ii] <- -opt$value
    my_val <- v2[ii]
  }

  ### Lower tail ...
  par_which <- mle_which
  my_val <- max_loglik
  ii <- 1
  sol <- mle[-which]
  while (my_val > conf_line){
    par_which <- par_which - inc
    opt <- stats::optim(sol, profiling_fn, ...)
    sol <- opt$par
    ii <- ii + 1
    x1[ii] <- par_which
    v1[ii] <- -opt$value
    my_val <- v1[ii]
  }

  # Find the limits of the confidence interval
  prof_lik <- c(rev(v1), v2)
  par_values <- c(rev(x1), x2)
  # Find where the curve crosses conf_line
  temp <- diff(prof_lik - conf_line > 0)
  # Find the upper limit of the confidence interval
  loc <- which(temp == -1)
  x1 <- par_values[loc]
  x2 <- par_values[loc + 1]
  y1 <- prof_lik[loc]
  y2 <- prof_lik[loc + 1]
  up_lim <- x1 + (conf_line - y1) * (x2 - x1) / (y2 - y1)
  # Find the lower limit of the confidence interval
  loc <- which(temp == 1)
  x1 <- par_values[loc]
  x2 <- par_values[loc+1]
  y1 <- prof_lik[loc]
  y2 <- prof_lik[loc+1]
  low_lim <- x1 + (conf_line - y1) * (x2 - x1) / (y2 - y1)
  par_prof <- c(lower = low_lim, mle_which, upper = up_lim)
  return(list(par_prof = par_prof, crit = conf_line,
              for_plot = cbind(par_values = par_values,
                               prof_loglik = prof_lik)))
}

# ======================== Faster GEV profiling function ==================== #

#' @keywords internal
#' @rdname profileCI-internal
faster_profile_ci <- function(negated_loglik_fn, which = 1, level, mle,
                              ci_sym_mat, inc, ...) {

  # The -log-likelihood to profile over parameters in par other than par[which]
  profiling_fn <- function(par, ...) {
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

  # Extract the sample maxima
  dots <- list(...)
  data <- dots$maxima_notNA$maxima

  ### Upper tail ...

  # We start from the upper limit of the symmetric confidence interval, using
  # gev_profile_init() to set initial estimates of the GEV parameters other
  # than the parameter which

  # Call gev_profile_init()
  if (which == 1) {
    par_which <- ci_sym_mat[1, 2]
    init <- gev_profile_init(data = data, mu = par_which)[-1]
  } else if (which == 2) {
    par_which <- ci_sym_mat[2, 2]
    init <- gev_profile_init(data = data, sigma = par_which)[-2]
  } else {
    par_which <- ci_sym_mat[3, 2]
    init <- gev_profile_init(data = data, xi = par_which)[-3]
  }
  # Calculate the profile log-likelihood at the initial values
  # If this is greater than conf_line then we search upwards
  # If this is smaller than conf_line then we search downwards

  ii <- 2
  opt <- stats::optim(init, profiling_fn, ...)
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
    opt <- stats::optim(sol, profiling_fn, ...)
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

  ### Lower tail ...

  # We start from the lower limit of the symmetric confidence interval, using
  # gev_profile_init() to set initial estimates of the GEV parameters other
  # than the parameter which

  # Call gev_profile_init()
  if (which == 1) {
    par_which <- ci_sym_mat[1, 1]
    init <- gev_profile_init(data = data, mu = par_which)[-1]
  } else if (which == 2) {
    par_which <- ci_sym_mat[2, 1]
    init <- gev_profile_init(data = data, sigma = par_which)[-2]
  } else {
    par_which <- ci_sym_mat[3, 1]
    init <- gev_profile_init(data = data, xi = par_which)[-3]
  }
  # Calculate the profile log-likelihood at the initial values
  # If this is greater than conf_line then we search downwards
  # If this is smaller than conf_line then we search upwards

  ii <- 2
  opt <- stats::optim(init, profiling_fn, ...)
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
    opt <- stats::optim(sol, profiling_fn, ...)
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

  # Find the limits of the confidence interval
  prof_lik <- c(rev(v1), v2)
  par_values <- c(rev(x1), x2)
  # Find where the curve crosses conf_line
  temp <- diff(prof_lik - conf_line > 0)
  # Find the upper limit of the confidence interval
  loc <- which(temp == -1)
  x1 <- par_values[loc]
  x2 <- par_values[loc + 1]
  y1 <- prof_lik[loc]
  y2 <- prof_lik[loc + 1]
  up_lim <- x1 + (conf_line - y1) * (x2 - x1) / (y2 - y1)
  # Find the lower limit of the confidence interval
  loc <- which(temp == 1)
  x1 <- par_values[loc]
  x2 <- par_values[loc+1]
  y1 <- prof_lik[loc]
  y2 <- prof_lik[loc+1]
  low_lim <- x1 + (conf_line - y1) * (x2 - x1) / (y2 - y1)
  par_prof <- c(lower = low_lim, mle_which, upper = up_lim)
  return(list(par_prof = par_prof, crit = conf_line,
              for_plot = cbind(par_values = par_values,
                               prof_loglik = prof_lik)))
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
