## True effects, random-number streams, replication runner, performance measures.

test_that("tau_super has the analytic values (spec 4)", {
  expect_equal(tau_super("sim1")$tau_super, c(-1.5, -1, -1.5, 0.5, 0, -0.5))
  expect_equal(tau_super("sim2")$tau_super, rep(0, 36))
})

test_that("Monte Carlo population moments agree with known E[W] and tau_super", {
  mom <- population_moments("sim1", n_chunks = 4L, chunk_size = 2e5, seed = 7L, cores = 2L)
  expect_equal(sum(mom$pi), 1)
  ## pi-weighted average of E[W | T = t] is E[W]
  expect_equal(unname(colSums(mom$EW_T * mom$pi)), unname(E_W), tolerance = 0.01)
  tr <- population_truth(mom)
  expect_equal(tr$tau_super, c(-1.5, -1, -1.5, 0.5, 0, -0.5))
  expect_true(all(tr$se_tau_mix > 0 & tr$se_tau_mix < 0.05))
  expect_true(all(tr$varsum > 0))
})

test_that("replication streams are reproducible and distinct", {
  s1 <- replication_seeds(123L, 3); s2 <- replication_seeds(123L, 3)
  expect_identical(s1, s2)
  expect_false(identical(s1[[1]], s1[[2]]))
  r_a <- run_replication(SIMS$sim1, 2, s1[[2]])
  r_b <- run_replication(SIMS$sim1, 2, s2[[2]])
  expect_identical(r_a$m, r_b$m)
  r_c <- run_replication(SIMS$sim1, 3, s1[[3]])
  expect_false(isTRUE(all.equal(r_a$m, r_c$m)))
  expect_equal(dim(r_a$m), c(length(ESTIMATORS), 4L))
})

test_that("run_simulation saves each replication and resumes", {
  out <- file.path(tempdir(), "resume_test"); unlink(out, recursive = TRUE)
  done1 <- run_simulation(SIMS$sim1, 1:2, out, cores = 2L, log_file = tempfile())
  expect_setequal(done1, 1:2)
  first <- readRDS(rep_file(out, 1))
  done2 <- run_simulation(SIMS$sim1, 1:3, out, cores = 2L, log_file = tempfile())
  expect_equal(done2, 3L)                                 # only the missing one ran
  expect_identical(readRDS(rep_file(out, 1)), first)
  ## a replication run alone equals the same replication inside a batch
  solo <- run_replication(SIMS$sim1, 3, replication_seeds(SIMS$sim1$master_seed, 3)[[3]])
  expect_identical(readRDS(rep_file(out, 3))$m, solo$m)
  col <- collect_simulation(out)
  expect_equal(nrow(col$estimates), 3 * length(ESTIMATORS) * 6)
})

test_that("performance measures match hand computation", {
  est <- data.frame(sim = "x", rep = 1:4, estimator = "DIF", t = 1, tprime = 2,
                    est = c(1, 2, 3, 6))
  truth <- data.frame(t = 1, tprime = 2, tau_super = 2)
  p <- performance(est, truth)
  expect_equal(p$bias, 1); expect_equal(p$abs_bias, 1)
  expect_equal(p$mse, mean(c(1, 0, 1, 16)))
  expect_equal(p$sd, sd(c(1, 2, 3, 6))); expect_equal(p$mcse_bias, sd(c(1, 2, 3, 6)) / 2)
})
