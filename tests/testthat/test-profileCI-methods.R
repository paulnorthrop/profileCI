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

# Repeat for the nls example

DNase1 <- subset(datasets::DNase, Run == 1)
fm1DNase1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)
prof1 <- profileCI(fm1DNase1, parm = 2)
prof2 <- profileCI(fm1DNase1, parm = "xmid")

test_that("Profile-based intervals for non-linear least squares regression", {
  expect_equal(prof1, prof2)
})
