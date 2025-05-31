# 1. Check that my symmetric intervals agree with confint.default

my_sym <- profileCI(glm.D93, loglik = poisson_loglik, profile = FALSE)
sym <- confint.default(glm.D93)

test_that("Symmetric intervals for Poisson GLM", {
  expect_equal(my_sym, sym, ignore_attr = TRUE)
})

