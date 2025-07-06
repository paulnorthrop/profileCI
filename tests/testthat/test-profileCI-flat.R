# 1. Check that if the profile log-likelihood doesn't drop to the cutoff then
#    +/- Inf is returned

# Use the data and example from
# https://github.com/paulnorthrop/profileCI/issues/1

x <- c(
  rep(-13.0, 6), rep(-12.0, 6), rep(-11.0, 6), rep(-10.0, 6), rep(-9.5, 6),
  rep(-9.0, 6), rep(-8.5, 6), rep(-8.0, 6), rep(-7.5, 6), rep(-7.0, 6),
  rep(-6.0, 6), rep(-5.0, 6)
)

w <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
       1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1)

y <- c(-0.035485178,-0.401627191,-0.094542019,-0.018607493,-0.185361946,
       0.072983324,-0.030225560,-0.218935722,-0.157890669,-0.049344803,
       -0.191411655,-0.166177041,0.010855507,-0.248252576,-0.115626745,
       -0.074683133,-0.166047692,-0.171485360,0.024284778,0.717596342,
       -0.180014600,0.027237054,-0.145052963,-0.105313167,-0.002500474,
       -0.184722343,-0.024193756,-0.003709391,-0.122355353,-0.105541012,
       -0.008803394,0.034549781,-0.121209258,-0.053846704,-0.165498737,
       -0.141773318,-0.007910983,0.677203125,-0.177188990,0.029961310,
       -0.164899368,-0.129556913,-0.035161841,0.205791107,-0.143891442,
       -0.053973286,-0.185524392,-0.131626915,-0.008794771,-0.270168575,
       -0.166917164,-0.157324396,-0.121862414,-0.131680240,-0.037675249,
       -0.565056503,-0.119697685,-0.101336812,-0.174113970,-0.123225713,
       -0.069556295,-2.046112059,-0.163627775,-0.081491021,-0.197287714,
       -0.153810355,-0.180387660,-2.514209890,-0.125615150,-0.140120308,
       -0.189019779,-0.170729106)

logistic4_fn <- function(x, theta) {
  alpha <- theta[1]
  delta <- theta[2]
  phi <- theta[3]
  eta <- 1
  return(alpha + delta / (1 + exp(-eta * (x - phi))))
}

ll_drda <- function(pars, x, y, w, ...) {
  # Calculate predicted values
  mu <- logistic4_fn(x, pars)
  # Calculate weighted residual sum of squares
  res <- sqrt(w) * (y - mu)
  rss <- sum(res^2)
  n <- length(y)
  # Log-likelihood based on RSS (matches drda's approach)
  ll <- -0.5 * n * (log(2 * pi * rss / n) + 1)
  return(ll)
}

coefs <- c(alpha = -0.09610808, delta = -0.08238220, phi = -6.66472838)

vc <- matrix(
  c(0.0003491587, -0.0002983981, -0.0132819893,
    -0.0002983981, 0.0037652802, -0.0606452437,
    -0.0132819893,-0.0606452437, 3.0571495500),
  dimnames = list(
    c("alpha", "delta", "phi"),
    c("alpha", "delta", "phi")
  ),
  nrow = 3, byrow = TRUE,
)

dummy2 <- list()
dummy2$coefficients <- coefs
dummy2$vcov <- vc
class(dummy2) <- "foo"
coef.foo <- function(x) x$coefficients
vcov.foo <- function(x) x$vcov

# 70% CI for phi
p70false <- profileCI(dummy2, loglik = ll_drda, parm = "phi", level = 0.7,
                      x = x, y = y, w = w, faster = FALSE)
p70true <- profileCI(dummy2, loglik = ll_drda, parm = "phi", level = 0.7,
                     x = x, y = y, w = w, faster = TRUE)
test_that("Flat upper profile log-likelihood, faster = FALSE", {
  expect_equal(p70false[2], Inf)
})
test_that("Flat upper profile log-likelihood, faster = TRUE", {
  expect_equal(p70true[2], Inf)
})

# 95% CI for phi
p95false <- profileCI(dummy2, loglik = ll_drda, parm = "phi", level = 0.95,
                      x = x, y = y, w = w, faster = FALSE)
p95true <- profileCI(dummy2, loglik = ll_drda, parm = "phi", level = 0.95,
                     x = x, y = y, w = w, faster = TRUE)
test_that("Flat lower and upper profile log-likelihood, faster = FALSE", {
  expect_equal(p95false, c(-Inf, Inf), ignore_attr = TRUE)
})
test_that("Flat lower and upper profile log-likelihood, faster = TRUE", {
  expect_equal(p95true, c(-Inf, Inf), ignore_attr = TRUE)
})

