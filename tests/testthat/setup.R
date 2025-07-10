# Create the GLM example on which tests are based

## From example(glm)
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3, 1, 9); treatment <- gl(3, 3)
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())
# Fit the intercept only model
glm.0 <- glm(counts ~ 1, family = poisson())

poisson_loglik <- function(par) {
  lambda <- exp(model.matrix(glm.D93) %*% par)
  loglik <- stats::dpois(x = glm.D93$y, lambda = lambda, log = TRUE)
  return(sum(loglik))
}

poisson_loglik_0 <- function(par) {
  lambda <- exp(model.matrix(glm.0) %*% par)
  loglik <- stats::dpois(x = glm.0$y, lambda = lambda, log = TRUE)
  return(sum(loglik))
}

## Create data used in the nls() example

Run <- rep_len(1, 16)
conc <- c(0.04882812, 0.04882812, 0.19531250, 0.19531250, 0.39062500,
          0.39062500, 0.78125000, 0.78125000, 1.56250000, 1.56250000,
          3.12500000, 3.12500000, 6.25000000, 6.25000000, 12.50000000,
          12.50000000)
density <- c(0.017, 0.018, 0.121, 0.124, 0.206, 0.215, 0.377, 0.374, 0.614,
             0.609, 1.019, 1.001, 1.334, 1.364, 1.730, 1.710)
DNase1 <- data.frame(Run = Run, conc = conc, density = density)
