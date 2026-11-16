## Multinomial logit by Newton-Raphson.

test_that("with two classes it equals binary logistic regression (glm)", {
  set.seed(1); n <- 300; X <- cbind(a = rnorm(n), b = runif(n))
  y <- rbinom(n, 1, plogis(-.5 + X %*% c(1, -2)))
  f <- fit_mlogit(X, y + 1L, 2)
  g <- glm(y ~ X, family = binomial)
  expect_true(f$converged)
  expect_equal(f$fitted[, 2], unname(fitted(g)), tolerance = 1e-8)
  expect_equal(f$loglik, as.numeric(logLik(g)), tolerance = 1e-10)
})

test_that("agrees with nnet::multinom and attains at least its log-likelihood", {
  skip_if_not_installed("nnet")
  set.seed(2); W <- draw_covariates(1500); T <- draw_treatment(true_gps(W, "sim2"), runif(1500))
  X <- gps_design(W, "main")
  f <- fit_mlogit(X, T, 9)
  g <- nnet::multinom(factor(T) ~ X, trace = FALSE, maxit = 2000, reltol = 1e-14)
  expect_true(f$converged)
  expect_equal(unname(f$fitted), unname(fitted(g)), tolerance = 1e-4)
  expect_gte(f$loglik, as.numeric(logLik(g)) - 1e-6)
})

test_that("score equations hold at the solution", {
  set.seed(3); W <- draw_covariates(800); T <- draw_treatment(true_gps(W, "sim1"), runif(800))
  X <- cbind(1, gps_design(W, "cs"))
  f <- fit_mlogit(X[, -1], T, 4)
  Yind <- outer(T, 1:4, "==") * 1
  expect_lt(max(abs(crossprod(X, Yind - f$fitted))), 1e-6)
})

test_that("recovers known coefficients in a large sample", {
  set.seed(4); n <- 40000; X <- cbind(rnorm(n), rnorm(n))
  B <- rbind(c(.5, -1), c(1, .5), c(-.5, 0))        # (intercept, x1, x2) for k = 2, 3
  P <- softmax_rows(cbind(0, cbind(1, X) %*% B))
  T <- draw_treatment(P, runif(n))
  f <- fit_mlogit(X, T, 3)
  ## back-transform the standardised coefficients
  b <- f$coef_std[-1, , drop = FALSE] / f$scale
  a <- f$coef_std[1, ] - colSums(b * f$center)
  expect_equal(unname(rbind(a, b)), B, tolerance = 0.05)
})

test_that("fitted probabilities do not depend on column scaling", {
  set.seed(5); W <- draw_covariates(600); T <- draw_treatment(true_gps(W, "sim1"), runif(600))
  X <- gps_design(W, "main")
  expect_equal(fit_mlogit(X, T, 4)$fitted, fit_mlogit(X %*% diag(c(1, 10, .1, 3, 1, 5, 2)) + 7, T, 4)$fitted,
               tolerance = 1e-8)
})

test_that("GPS design matrices have the specified columns", {
  W <- draw_covariates(5)
  expect_equal(colnames(gps_design(W, "main")), COV_NAMES)
  expect_equal(colnames(gps_design(W, "cs")),
               c(COV_NAMES, "x1_sq", "x2_sq", "x3_sq", "x4_sq", "z1_sq", "z2_sq"))
})
