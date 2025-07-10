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

# 5. Tests for an intercept only Poisson GLM

conf0 <- confint(glm.0, parm = 1)
prof0 <- profileCI(glm.0, loglik = poisson_loglik_0, profile = TRUE,
                   mult = 32, faster = TRUE)
test_that("Profile-based intervals for intercept only Poisson GLM", {
  expect_equal(conf0, prof0, tolerance = 1e-5, ignore_attr = TRUE)
})

# 6. nls() fast = TRUE vs fast = FALSE

# From example(nls)
test_that("nls() example, fast vs slow", {
  Run <- rep_len(1, 16)
  conc <- c(0.04882812, 0.04882812, 0.19531250, 0.19531250, 0.39062500,
            0.39062500, 0.78125000, 0.78125000, 1.56250000, 1.56250000,
            3.12500000, 3.12500000, 6.25000000, 6.25000000, 12.50000000,
            12.50000000)
  density <- c(0.017, 0.018, 0.121, 0.124, 0.206, 0.215, 0.377, 0.374, 0.614,
               0.609, 1.019, 1.001, 1.334, 1.364, 1.730, 1.710)
  DNase1 <- data.frame(Run = Run, conc = conc, density = density)
  fm1DNase1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal),
                   data = DNase1)
  prof1 <- profileCI(fm1DNase1, faster = TRUE)
  prof2 <- profileCI(fm1DNase1, faster = FALSE)
  expect_equal(prof1, prof2, tolerance = 1e-5, ignore_attr = TRUE)
})
