# Create the GLM example on which tests are based

## From example(glm)
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3, 1, 9); treatment <- gl(3, 3)
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())
# Fir the intercept only model
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
