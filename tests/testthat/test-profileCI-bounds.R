# Check that imposing lower and upper bounds lb and ub results in the correct
# confidence limits becoming NA

## faster = FALSE

# Bounded search
prof <- profileCI(glm.D93, parm = 1:4, faster = FALSE,
                  lb = rep(-0.6, 4), ub = rep(3.2, 4))
# Elements [2, 1], [3, 1] and [1, 2] of prof should become NA in prof_bounded
theNAs <- c(prof[2, 1], prof[3, 1], prof[1, 2])

test_that("Bounds lb and ub work, faster = FALSE", {
  expect_true(all(is.na(theNAs)))
})

## faster = TRUE

# Bounded search
prof <- profileCI(glm.D93, parm = names(coef(glm.D93))[1:4], faster = TRUE,
                  lb = rep(-0.6, 4), ub = rep(3.2, 4))
# Elements [2, 1], [3, 1] and [1, 2] of prof should become NA in prof_bounded
theNAs <- c(prof[2, 1], prof[3, 1], prof[1, 2])

test_that("Bounds lb and ub work, faster = TRUE", {
  expect_true(all(is.na(theNAs)))
})
