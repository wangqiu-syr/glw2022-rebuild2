## CM metric variants, matching details and the Stage-4 variance estimators.

toy <- list(W = cbind(x = c(0, 2, 4, 1, 3, 5), z = c(0, 0, 1, 0, 1, 1)),
            T = c(1L, 1L, 1L, 2L, 2L, 2L), Y = c(1, 3, 5, 10, 20, 30))

test_that("stacked metric equals one Mahalanobis distance on all covariates", {
  set.seed(1); W <- draw_covariates(60)
  B <- cm_blocks(W, "stacked", COV_NAMES[1:4], COV_NAMES[5:7])
  V <- crossprod(sweep(W, 2, colMeans(W))) / nrow(W)
  D <- cm_dist(B, 1:10)
  ref <- sapply(1:10, function(j) sqrt(stats::mahalanobis(W, W[j, ], V)))
  expect_equal(unname(D), unname(ref), tolerance = 1e-10)
})

test_that("x_mahal_z_euclid adds a Euclidean Z distance to the X Mahalanobis distance", {
  set.seed(2); W <- draw_covariates(40)
  Vx <- crossprod(sweep(W[, 1:4], 2, colMeans(W[, 1:4]))) / 40
  D <- cm_dist(cm_blocks(W, "x_mahal_z_euclid", COV_NAMES[1:4], COV_NAMES[5:7]), 1:5)
  ref <- sapply(1:5, function(j) sqrt(stats::mahalanobis(W[, 1:4], W[j, 1:4], Vx)) +
                  sqrt(rowSums(sweep(W[, 5:7], 2, W[j, 5:7])^2)))
  expect_equal(unname(D), unname(ref), tolerance = 1e-10)
})

test_that("metric choice can change the match", {
  ## x and z strongly correlated: the stacked metric discounts moving along the
  ## correlation, the block metric does not.
  W <- cbind(x = c(0, 1, 2, 3, 0.9, 2.1), z = c(0, 1, 2, 3, 0, 2.1))
  T <- c(1L, 1L, 1L, 1L, 2L, 2L); Y <- c(1, 2, 3, 4, 0, 0)
  blk <- est_cm(W, T, Y, 2, x_cols = "x", z_cols = "z", metric = "blocks")$Yhat[5, 1]
  stk <- est_cm(W, T, Y, 2, x_cols = "x", z_cols = "z", metric = "stacked")$Yhat[5, 1]
  expect_false(blk == stk)
  expect_equal(est_cm(W, T, Y, 2, x_cols = "x", z_cols = "z")$Yhat,
               est_cm(W, T, Y, 2, x_cols = "x", z_cols = "z", metric = "blocks")$Yhat)
})

test_that("match details reproduce the hand-computed toy example", {
  ## Level 1 (donors 1-3): matches 1->1, 2->2, 3->3, 4->{1,2}, 5->3, 6->3, so
  ## w = (1.5, 1.5, 3). Within-group nearest donors: 1<->2, 3->2, so
  ## sig2 = ((1-3)^2/2, (3-1)^2/2, (5-3)^2/2) = (2, 2, 2).
  ## Level 2 (donors 4-6): matches 1->4, 2->4, 3->{5,6}, so w = (3, 1.5, 1.5);
  ## nearest donors 4->5, 5<->6, so sig2 = (50, 50, 50).
  r <- est_cm(toy$W, toy$T, toy$Y, 2, x_cols = "x", z_cols = "z", details = TRUE)
  d1 <- r$details[[1]]; d2 <- r$details[[2]]
  expect_equal(d1$w, c(1.5, 1.5, 3)); expect_equal(d1$sig2, c(2, 2, 2))
  expect_equal(d1$v, c(2, 2, 2, 1, 2, 2))              # unit 4: .25*2 + .25*2
  expect_equal(d2$w, c(3, 1.5, 1.5)); expect_equal(d2$sig2, c(50, 50, 50))
  expect_equal(d2$v, c(50, 50, 25, 50, 50, 50))        # unit 3: .25*50 + .25*50
})

test_that("pair variances reproduce the hand-computed toy example", {
  ## var_noise = (sum w^2 sig2 over both groups) / n^2 = (27 + 675) / 36 = 19.5
  ## unit effects tau_i = (-9, -7, -20, -8, -15, -25), noise v_i = (52, 52, 27, 51, 52, 52)
  ## group 1: var 49 - mean v 43.667 = 5.333; group 2: var 73 - mean v 51.667 = 21.333
  ## var_sate = (3 * 5.333 + 3 * 21.333) / 36 = 80 / 36
  r <- est_cm(toy$W, toy$T, toy$Y, 2, x_cols = "x", z_cols = "z", details = TRUE)
  pv <- pair_variances_matching(r, toy$T, 2)
  expect_equal(pv$var_noise, 19.5)
  expect_equal(pv$var_sate, 80 / 36)
})

test_that("DIF variance is s_t^2/n_t + s_t'^2/n_t'", {
  pv <- pair_variances_dif(toy$Y, toy$T, 2)
  expect_equal(pv$var_full, var(c(1, 3, 5)) / 3 + var(c(10, 20, 30)) / 3)
})

test_that("GPS details use within-group nearest neighbours on the score", {
  S <- cbind(c(.6, .3, .55, .2), c(.4, .7, .45, .8)); T <- c(1L, 1L, 2L, 2L); Y <- c(1, 2, 10, 20)
  r <- est_gps(S, T, Y, 2, details = TRUE)
  expect_equal(r$details[[1]]$sig2, c(.5, .5))         # (1 - 2)^2 / 2
  expect_equal(r$details[[1]]$w, c(2, 2))              # units 1, 3 -> 1; units 2, 4 -> 2
  expect_equal(r$details[[2]]$sig2, c(50, 50))
})

test_that("variance estimator is calibrated in a simple randomised setting", {
  ## Two groups of 150, no confounding, heterogeneous effect: tau = 1.
  sim_once <- function() {
    n <- 300; T <- rep(1:2, each = 150)
    W <- cbind(matrix(rnorm(4 * n), n), matrix(rnorm(3 * n), n)); colnames(W) <- COV_NAMES
    Y <- ifelse(T == 1, 1 + W[, 1] + W[, 5], 2 + 2 * W[, 1] + W[, 5]) + rnorm(n)
    cm <- est_cm(W, T, Y, 2, details = TRUE); S <- cbind(plogis(W[, 1]), plogis(W[, 1]))
    gp <- est_gps(S, T, Y, 2, details = TRUE)
    pc <- pair_variances_matching(cm, T, 2); pg <- pair_variances_matching(gp, T, 2)
    c(cm = diff(rev(cm$m)), cm_v = pc$var_noise + pc$var_sate,
      gps = diff(rev(gp$m)), gps_v = pg$var_noise + pg$var_sate)
  }
  set.seed(9); R <- t(replicate(200, sim_once()))
  expect_gt(mean(R[, "cm_v"]) / var(R[, "cm"]), 0.7); expect_lt(mean(R[, "cm_v"]) / var(R[, "cm"]), 1.45)
  expect_gt(mean(R[, "gps_v"]) / var(R[, "gps"]), 0.7); expect_lt(mean(R[, "gps_v"]) / var(R[, "gps"]), 1.45)
})

test_that("extra runs reproduce the Stage-3 data sets and estimates", {
  cfg <- SIMS$sim1; s <- replication_seeds(cfg$master_seed, 2)[[2]]
  base <- run_replication(cfg, 2, s)
  ext <- run_replication_cm_metrics(cfg, 2, s)
  expect_equal(ext$m["CM", ], base$m["CM", ], tolerance = 1e-12)
  s4 <- run_replication_stage4(cfg, 2, s)
  expect_equal(s4$m, base$m, tolerance = 1e-12)
  expect_true(all(s4$var$var_full > 0))
  expect_equal(nrow(s4$var), length(ESTIMATORS) * 6)
})
