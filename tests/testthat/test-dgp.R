## Data-generating process: hand-computed values and distributional checks.

w0 <- matrix(c(2, 2, -1, 0.5, -2, 3, 1), 1, dimnames = list(NULL, COV_NAMES))

test_that("eta_t matches hand computation at a fixed point (Sim 1 and Sim 2)", {
  ## w0 = (x1..x4, z1..z3) = (2, 2, -1, 0.5, -2, 3, 1)
  ## eta2 = .1(4 + 2 - 1 + .5 - 2 + 3 + 1) = .75
  ## eta3 = .1(2 + 4 + 1 + .5 - 2 + 3 + 1) = .95
  ## eta4 = .1(2 + 4 + 1 + .25 - 2 + 3 + 1) = .925
  expect_equal(as.vector(eta_matrix(w0, "sim1")), c(0, .75, .95, .925))
  ## eta5 = .1(2+4+1+.25+4+3+1) = 1.525; eta6 = .1(2+4+1+.25+4+9+1) = 2.125
  ## eta7 = eta6 (z3^2 = z3 = 1); eta8 = .1(4+4+1+.25+4+9+1) = 2.325
  ## eta9 = .1(2+2-1+.5-2+3+1) = .55
  expect_equal(as.vector(eta_matrix(w0, "sim2")),
               c(0, .75, .95, .925, 1.525, 2.125, 2.125, 2.325, .55))
})

test_that("true GPS is the softmax of eta and rows sum to one", {
  e <- c(0, .75, .95, .925)
  expect_equal(as.vector(true_gps(w0, "sim1")), exp(e) / sum(exp(e)))
  W <- draw_covariates(500)
  for (s in c("sim1", "sim2")) expect_equal(rowSums(true_gps(W, s)), rep(1, 500))
  ## extreme linear predictors do not overflow
  expect_equal(as.vector(softmax_rows(matrix(c(0, 800, 799), 1))),
               c(0, 1, exp(-1)) / (1 + exp(-1)), tolerance = 1e-12)
})

test_that("outcome means mu_t match hand computation", {
  ## Sim 1, t = 2: 2*2 + 3*2 + 1*(-1) + 2*.5 + 2*(-2) + 2*3 + 2*1 = 14
  expect_equal(mu_matrix(w0, "sim1")[1, 2], 14)
  ## Sim 2, t = 1: 2 + 2 - 1 + .5 - 2 + (3 - 1) + (1 - .5) = 4
  expect_equal(mu_matrix(w0, "sim2")[1, 1], 4)
  ## Sim 2, t = 6: 6*2 + 3*2 - 1 + 2*.5 + 2*(-2) + 2*(3 - 1) + 9*(1 - .5) = 22.5
  expect_equal(mu_matrix(w0, "sim2")[1, 6], 22.5)
  ## Sim 2, t = 8: 8*2 + 2 + 5*(-1) + 2.5*.5 + 4*(-2) + (3 - 1) + 7*.5 = 11.75
  expect_equal(mu_matrix(w0, "sim2")[1, 8], 11.75)
})

test_that("draw_treatment inverts the cumulative probabilities", {
  P <- matrix(c(.2, .3, .5), 4, 3, byrow = TRUE)
  expect_equal(draw_treatment(P, c(.1, .25, .5, .51)), c(1L, 2L, 2L, 3L))
  expect_equal(draw_treatment(matrix(c(.4, .6), 2, 2, byrow = TRUE), c(.4, .41)), c(1L, 2L))
})

test_that("covariate distribution has the specified moments", {
  set.seed(11); N <- 4e5; W <- draw_covariates(N)
  expect_equal(unname(cov(W[, 1:4])), SIGMA_X, tolerance = 0.02)
  expect_equal(unname(colMeans(W)), unname(E_W), tolerance = 0.01)
  expect_true(all(W[, "z1"] >= -3 & W[, "z1"] <= 3))
  expect_equal(var(W[, "z1"]), 3, tolerance = 0.02)          # Uniform(-3,3): 36/12
  expect_equal(var(W[, "z2"]), 2, tolerance = 0.03)          # chi2_1
  expect_setequal(unique(W[, "z3"]), c(0, 1))
  expect_equal(unname(cor(W[, 1:4], W[, 5:7])), matrix(0, 4, 3), tolerance = 0.01)
})

test_that("simulate_sample gives exactly n_t per group and Y = mu_T + eps", {
  for (s in names(SIMS)) {
    cfg <- SIMS[[s]]; set.seed(3); d <- simulate_sample(cfg)
    expect_equal(as.vector(table(d$T)), rep(cfg$n_per_group, cfg$A))
    expect_equal(nrow(d$W), cfg$A * cfg$n_per_group)
    res <- d$Y - d$mu[cbind(seq_along(d$T), d$T)]
    expect_lt(abs(mean(res)), 0.15); expect_lt(abs(sd(res) - 1), 0.1)
    expect_equal(d$P_true, true_gps(d$W, s))
  }
})

test_that("the group-filling scheme samples from W | T = t", {
  ## Compare group means of a large filled sample with the weighting identity
  ## E[W | T = t] = E[W p_t(W)] / E[p_t(W)] computed from independent draws.
  cfg <- modifyList(SIMS$sim1, list(n_per_group = 20000L))
  set.seed(5); d <- simulate_sample(cfg, batch_size = 5000L)
  set.seed(6); W <- draw_covariates(1e6); P <- true_gps(W, "sim1")
  for (t in 1:4) {
    target <- colSums(W * P[, t]) / sum(P[, t])
    got <- colMeans(d$W[d$T == t, ])
    se <- apply(d$W[d$T == t, ], 2, sd) / sqrt(cfg$n_per_group)
    expect_true(all(abs(got - target) < 4.5 * se), info = paste("t =", t))
  }
})

test_that("Sim 2 treatments 6 and 7 have identical propensity functions (z3^2 = z3)", {
  W <- draw_covariates(1000)
  P <- true_gps(W, "sim2")
  expect_equal(P[, 6], P[, 7])
})
