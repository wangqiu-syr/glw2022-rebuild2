## Performance measures (docs/spec.md section 5).
## est: long data frame (sim, rep, estimator, t, tprime, est);
## truth: data frame (t, tprime, <truth_col>).
performance <- function(est, truth, truth_col = "tau_super") {
  d <- merge(est, truth[, c("t", "tprime", truth_col)], by = c("t", "tprime"))
  d$err <- d$est - d[[truth_col]]
  key <- interaction(d$estimator, d$t, d$tprime, drop = TRUE)
  out <- do.call(rbind, lapply(split(d, key), function(g) {
    R <- nrow(g)
    data.frame(sim = g$sim[1], estimator = g$estimator[1], t = g$t[1], tprime = g$tprime[1],
               truth = truth_col, tau = g[[truth_col]][1], R = R,
               mean_est = mean(g$est), bias = mean(g$err), abs_bias = abs(mean(g$err)),
               mse = mean(g$err^2), sd = sd(g$est),
               mcse_bias = sd(g$est) / sqrt(R), mcse_mse = sd(g$err^2) / sqrt(R))
  }))
  out$estimator <- factor(out$estimator, levels = c(ESTIMATORS, setdiff(unique(out$estimator), ESTIMATORS)))
  out <- out[order(out$estimator, out$t, out$tprime), ]
  out$estimator <- as.character(out$estimator)
  rownames(out) <- NULL
  out
}
