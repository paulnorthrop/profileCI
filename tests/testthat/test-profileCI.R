# 1. Check that my symmetric intervals agree with confint.default()

my_sym <- profileCI(glm.D93, loglik = poisson_loglik, profile = FALSE)
sym <- confint.default(glm.D93)

test_that("Symmetric intervals for Poisson GLM", {
  expect_equal(my_sym, sym, ignore_attr = TRUE)
})

# 2. Check that profile-based intervals are close enough to confint.glm()

prof <- confint(glm.D93, parm = 1)

# faster = TRUE
my_prof_faster <- profileCI(glm.D93, loglik = poisson_loglik, parm = 1,
                            profile = TRUE, mult = 32, faster = TRUE)
test_that("Profile-based intervals for Poisson GLM, faster = TRUE", {
  expect_equal(my_prof_faster, prof, tolerance = 1e-5, ignore_attr = TRUE)
})

# faster = FALSE
my_prof <- profileCI(glm.D93, loglik = poisson_loglik, parm = 1,
                     profile = TRUE, mult = 32, faster = FALSE)
test_that("Profile-based intervals for Poisson GLM, faster = FALSE", {
  expect_equal(my_prof, prof, tolerance = 1e-5, ignore_attr = TRUE)
})

# 3. Check that passing loglik and using logLikFn give the same results

my_prof_logLikFn <- profileCI(glm.D93, parm = 1,
                              profile = TRUE, mult = 32, faster = TRUE)
test_that("Profile-based intervals for Poisson GLM, faster = FALSE", {
  expect_equal(my_prof_faster, my_prof_logLikFn, tolerance = 1e-5,
               ignore_attr = FALSE)
})

# 4. Tests with epsilon > 0

my_prof_fast <- profileCI(glm.D93, loglik = poisson_loglik, parm = 2,
                          profile = TRUE, mult = 32, faster = TRUE,
                          epsilon = 1e-4)
my_prof_slow <- profileCI(glm.D93, loglik = poisson_loglik, parm = 2,
                          profile = TRUE, mult = 32, faster = FALSE,
                          epsilon = 1e-4)
test_that("Profile-based intervals for Poisson GLM, faster = TRUE, epsilon > 0", {
  expect_equal(my_prof_fast, my_prof_slow, tolerance = 1e-3, ignore_attr = TRUE)
})
