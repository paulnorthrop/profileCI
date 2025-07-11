# Check logLikFn()

## Create data used in the nls() example

Run <- rep_len(1, 16)
conc <- c(0.04882812, 0.04882812, 0.19531250, 0.19531250, 0.39062500,
          0.39062500, 0.78125000, 0.78125000, 1.56250000, 1.56250000,
          3.12500000, 3.12500000, 6.25000000, 6.25000000, 12.50000000,
          12.50000000)
density <- c(0.017, 0.018, 0.121, 0.124, 0.206, 0.215, 0.377, 0.374, 0.614,
             0.609, 1.019, 1.001, 1.334, 1.364, 1.730, 1.710)
DNase1 <- data.frame(Run = Run, conc = conc, density = density)

fm1DNase1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)
logLikFn(fm1DNase1, pars = coef(fm1DNase1))

test_that("logLikFn.nls and logLik.nls agree at the MLE", {
  expect_equal(logLik(fm1DNase1), logLikFn(fm1DNase1, pars = coef(fm1DNase1)),
               ignore_attr = TRUE)
})
