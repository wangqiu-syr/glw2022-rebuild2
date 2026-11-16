## Estimators on small hand-made examples.

## Toy data: 6 units, A = 2, one X covariate and one Z covariate.
##   unit  T  x  z   Y
##     1   1  0  0   1
##     2   1  2  0   3
##     3   1  4  1   5
##     4   2  1  0  10
##     5   2  3  1  20
##     6   2  5  1  30
toy <- list(W = cbind(x = c(0, 2, 4, 1, 3, 5), z = c(0, 0, 1, 0, 1, 1)),
            T = c(1L, 1L, 1L, 2L, 2L, 2L), Y = c(1, 3, 5, 10, 20, 30))

test_that("DIF is the difference of group means", {
  m <- est_dif(toy$Y, toy$T, 2)
  expect_equal(m, c(3, 20))
  expect_equal(pairs_from_m(m)$est, -17)
})

test_that("Mahalanobis whitening reproduces (a'V^-1 a)^{1/2} with divisor n", {
  set.seed(1); M <- matrix(rnorm(40), 20, 2) %*% matrix(c(2, .5, 0, 1), 2)
  V <- crossprod(sweep(M, 2, colMeans(M))) / nrow(M)
  B <- mahalanobis_whiten(M)
  a <- M[3, ] - M[7, ]
  expect_equal(sqrt(sum((B[3, ] - B[7, ])^2)), sqrt(as.numeric(t(a) %*% solve(V) %*% a)))
  expect_equal(pairwise_dist(B, 7)[3, 1], sqrt(as.numeric(t(a) %*% solve(V) %*% a)))
  expect_identical(pairwise_dist(B, 1:20)[cbind(1:20, 1:20)], rep(0, 20))
})

test_that("impute_nn takes the nearest donor and averages exact ties", {
  D <- rbind(c(1, 2, 3),
             c(2, 2, 5),       # tie between donors 1 and 2
             c(4, 0, 0))       # tie between donors 2 and 3
  expect_equal(impute_nn(D, c(10, 20, 30)), c(10, 15, 25))
})

test_that("CM reproduces the hand-computed toy example", {
  ## sd_n(x) = sqrt(17.5/6), sd_n(z) = 0.5; distance = |dx|/sd_x + |dz|/0.5.
  ## Imputed Y(1): units 1-3 own; unit 4 ties units 1 and 2 -> 2; units 5, 6 -> unit 3 -> 5
  ## Imputed Y(2): unit 1 -> 4 (10); unit 2 -> 4 (10); unit 3 ties 5 and 6 -> 25; 4-6 own
  r <- est_cm(toy$W, toy$T, toy$Y, 2, x_cols = "x", z_cols = "z")
  expect_equal(r$Yhat[, 1], c(1, 3, 5, 2, 5, 5))
  expect_equal(r$Yhat[, 2], c(10, 10, 25, 10, 20, 30))
  expect_equal(r$m, c(21 / 6, 105 / 6))
  expect_equal(pairs_from_m(r$m)$est, -14)
})

test_that("CM uses the unit's own outcome for its own treatment", {
  set.seed(2); cfg <- SIMS$sim1; d <- simulate_sample(cfg)
  r <- est_cm(d$W, d$T, d$Y, cfg$A)
  expect_identical(r$Yhat[cbind(seq_along(d$T), d$T)], d$Y)
})

test_that("CM is invariant to affine transformations within each block", {
  set.seed(4); cfg <- SIMS$sim1; d <- simulate_sample(cfg)
  r1 <- est_cm(d$W, d$T, d$Y, cfg$A)
  W2 <- d$W
  W2[, 1:4] <- d$W[, 1:4] %*% matrix(c(1, .3, 0, 0, 0, 2, 0, 0, .5, 0, 1, 0, 0, 0, -1, 3), 4) + 5
  W2[, 5:7] <- d$W[, 5:7] %*% diag(c(2, .1, 7)) - 1
  r2 <- est_cm(W2, d$T, d$Y, cfg$A)
  expect_equal(r1$Yhat, r2$Yhat)
})

test_that("CM distance adds the X-block and Z-block distances", {
  ## X block alone would match unit 3 (x = 1.1) to unit 1 (x = 1);
  ## the Z block (z far apart) moves the match to unit 2.
  W <- cbind(x = c(1, 1.5, 1.1, 10), z = c(0, 5, 5, 5))
  T <- c(1L, 1L, 2L, 2L); Y <- c(100, 200, 0, 0)
  r <- est_cm(W, T, Y, 2, x_cols = "x", z_cols = "z")
  expect_equal(r$Yhat[3, 1], 200)
})

test_that("GPS matching reproduces a hand-computed example", {
  ##   unit T  s1   s2   Y
  ##   1    1  .60  .40  1
  ##   2    1  .30  .70  2
  ##   3    2  .55  .45  10
  ##   4    2  .20  .80  20
  S <- cbind(c(.6, .3, .55, .2), c(.4, .7, .45, .8))
  T <- c(1L, 1L, 2L, 2L); Y <- c(1, 2, 10, 20)
  ## Y(1) for unit 3: |.55 - .6| = .05 < |.55 - .3| = .25 -> 1
  ## Y(1) for unit 4: |.2 - .6| = .4 > |.2 - .3| = .1 -> 2
  ## Y(2) for unit 1: score s2 = .4; |.4 - .45| = .05, |.4 - .8| = .4 -> 10
  ## Y(2) for unit 2: s2 = .7; |.7 - .45| = .25, |.7 - .8| = .1 -> 20
  r <- est_gps(S, T, Y, 2)
  expect_equal(r$Yhat, cbind(c(1, 2, 1, 2), c(10, 20, 10, 20)))
  expect_equal(pairs_from_m(r$m)$est, 1.5 - 15)
})

test_that("GPS matches on the score of the target level only", {
  ## Changing scores for level 2 must not change the imputation of Y(1).
  set.seed(8); n <- 40; T <- rep(1:3, length.out = n); Y <- rnorm(n)
  S <- matrix(runif(3 * n), n); S2 <- S; S2[, 2:3] <- runif(2 * n)
  expect_equal(est_gps(S, T, Y, 3)$Yhat[, 1], est_gps(S2, T, Y, 3)$Yhat[, 1])
  expect_identical(est_gps(S, T, Y, 3)$Yhat[cbind(1:n, T)], Y)
})

test_that("pairwise estimates are differences of per-treatment means", {
  p <- pairs_from_m(c(1, 4, 9, 16))
  expect_equal(p$t, c(1, 1, 1, 2, 2, 3)); expect_equal(p$tprime, c(2, 3, 4, 3, 4, 4))
  expect_equal(p$est, c(-3, -8, -15, -5, -12, -7))
})
