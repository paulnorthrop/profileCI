# Check plot.profileCI()

my_prof <- profileCI(glm.D93, loglik = poisson_loglik, profile = TRUE,
                     mult = 32, faster = TRUE)

interval1 <- plot(my_prof, parm = 1)
interval2 <- plot(my_prof, parm = "outcome2")

test_that("Profile-based intervals for Poisson GLM", {
  expect_equal(interval1, my_prof[1, ])
})
test_that("Profile-based intervals for Poisson GLM", {
  expect_equal(interval2, my_prof["outcome2", ])
})
