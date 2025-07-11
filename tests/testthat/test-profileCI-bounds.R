# Check that imposing lower and upper bounds lb and ub results in the correct
# confidence limits becoming NA

## faster = FALSE

# Unbounded search
prof <- profileCI(glm.D93, faster = FALSE)
# Bounded search
prof_bounded <- profileCI(glm.D93, faster = FALSE,
                          lb = rep(-0.6, 5), ub = rep(3.2, 5))
# Elements [2, 1], [3, 1] and [1, 2] of prof should become NA in prof_bounded
prof[2, 1] <- NA
prof[3, 1] <- NA
prof[1, 2] <- NA

test_that("Bounds lb and ub work, faster = FALSE", {
  expect_equal(prof, prof_bounded, ignore_attr = TRUE)
})

## faster = TRUE

# Unbounded search
prof <- profileCI(glm.D93, faster = TRUE)
# Bounded search
prof_bounded <- profileCI(glm.D93, faster = TRUE,
                          lb = rep(-0.6, 5), ub = rep(3.2, 5))
# Elements [2, 1], [3, 1] and [1, 2] of prof should become NA in prof_bounded
prof[2, 1] <- NA
prof[3, 1] <- NA
prof[1, 2] <- NA

test_that("Bounds lb and ub work, faster = TRUE", {
  expect_equal(prof, prof_bounded, ignore_attr = TRUE)
})

