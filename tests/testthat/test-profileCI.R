# 1. Check that my symmetric intervals agree with confint.default()

my_sym <- profileCI(glm.D93, loglik = poisson_loglik, profile = FALSE)
sym <- confint.default(glm.D93)

test_that("Symmetric intervals for Poisson GLM", {
  expect_equal(my_sym, sym, ignore_attr = TRUE)
})

# 2. Check that profile-based intervals are close enough to confint.glm()

prof <- confint(glm.D93)

# faster = TRUE
my_prof <- profileCI(glm.D93, loglik = poisson_loglik, profile = TRUE,
                     mult = 32, faster = TRUE)
test_that("Profile-based intervals for Poisson GLM", {
  expect_equal(my_prof, prof, tolerance = 1e-5, ignore_attr = TRUE)
})

# faster = FALSE
my_prof <- profileCI(glm.D93, loglik = poisson_loglik, profile = TRUE,
                     mult = 32, faster = FALSE)
test_that("Profile-based intervals for Poisson GLM", {
  expect_equal(my_prof, prof, tolerance = 1e-5, ignore_attr = TRUE)
})
