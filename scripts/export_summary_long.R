## Exports rebuild 2's results in one tidy table for later comparisons
## (docs/comparison_interface.md). Reads only this repository's tracked result
## tables. From the repository root:  Rscript scripts/export_summary_long.R
## Writes results/rebuild2_summary_long.csv with columns
##   sim, estimator, t, tprime, measure, truth, tau, value, mcse, reps
EST_ORDER <- c("DIF", "CM", "GPS", "GPS_cs", "GPS_true", "CM_stacked", "CM_xmahal_zeuclid")
rows <- list()
add <- function(d) rows[[length(rows) + 1]] <<- d

for (s in c("sim1", "sim2")) {
  ## bias, |bias|, MSE, mean estimate, SD (main estimators and CM metric variants)
  pf <- read.csv(file.path("results", s, "performance.csv"))
  cmv <- read.csv(file.path("results/cm_metric_sensitivity", sprintf("performance_%s.csv", s)))
  pf <- rbind(pf, cmv[cmv$estimator != "CM", ])            # CM itself is already in performance.csv
  for (m in c("mean_est", "bias", "abs_bias", "mse")) {
    se <- switch(m, mse = pf$mcse_mse, pf$mcse_bias)
    add(data.frame(sim = s, estimator = pf$estimator, t = pf$t, tprime = pf$tprime, measure = m,
                   truth = if (m == "mean_est") NA else pf$truth, tau = if (m == "mean_est") NA else pf$tau,
                   value = pf[[m]], mcse = se, reps = pf$R))
  }
  one <- pf[pf$truth == "tau_super", ]                     # SD does not depend on the truth
  add(data.frame(sim = s, estimator = one$estimator, t = one$t, tprime = one$tprime, measure = "sd",
                 truth = NA, tau = NA, value = one$sd, mcse = one$sd / sqrt(2 * (one$R - 1)), reps = one$R))
  ## coverage and SE calibration (Stage 4)
  cv <- read.csv(file.path("results", s, "coverage.csv"))
  tr <- read.csv(sprintf("results/truth/%s_truth.csv", s))
  iv <- read.csv(gzfile(file.path("results", s, "intervals.csv.gz")))
  R <- length(unique(iv$rep))
  cv <- merge(cv, tr[, c("t", "tprime", "tau_super", "tau_mix")], by = c("t", "tprime"))
  pse <- function(p) sqrt(p * (1 - p) / R)
  add(data.frame(sim = s, estimator = cv$estimator, t = cv$t, tprime = cv$tprime, measure = "coverage",
                 truth = "tau_super", tau = cv$tau_super, value = cv$cov_super, mcse = pse(cv$cov_super), reps = R))
  add(data.frame(sim = s, estimator = cv$estimator, t = cv$t, tprime = cv$tprime, measure = "coverage",
                 truth = "tau_mix", tau = cv$tau_mix, value = cv$cov_mix, mcse = pse(cv$cov_mix), reps = R))
  ok <- !is.na(cv$cov_sate_cond)
  add(data.frame(sim = s, estimator = cv$estimator[ok], t = cv$t[ok], tprime = cv$tprime[ok],
                 measure = "coverage_cond", truth = "sample_effect", tau = NA, value = cv$cov_sate_cond[ok],
                 mcse = pse(cv$cov_sate_cond[ok]), reps = R))
  add(data.frame(sim = s, estimator = cv$estimator, t = cv$t, tprime = cv$tprime, measure = "coverage_oracle",
                 truth = "tau_super", tau = cv$tau_super, value = cv$cov_oracle_super,
                 mcse = pse(cv$cov_oracle_super), reps = R))
  iv$se <- sqrt(iv$var_full)
  ms <- aggregate(se ~ estimator + t + tprime, iv, function(x) c(m = mean(x), s = sd(x) / sqrt(length(x))))
  ms <- data.frame(ms[, 1:3], mean_se = ms$se[, "m"], mcse_se = ms$se[, "s"])
  add(data.frame(sim = s, estimator = ms$estimator, t = ms$t, tprime = ms$tprime, measure = "mean_se",
                 truth = NA, tau = NA, value = ms$mean_se, mcse = ms$mcse_se, reps = R))
  add(data.frame(sim = s, estimator = cv$estimator, t = cv$t, tprime = cv$tprime, measure = "se_ratio",
                 truth = NA, tau = NA, value = cv$se_ratio, mcse = NA, reps = R))
}
out <- do.call(rbind, rows)
out <- out[!duplicated(out[, c("sim", "estimator", "t", "tprime", "measure", "truth")]), ]
MEAS <- c("mean_est", "bias", "abs_bias", "mse", "sd", "coverage", "coverage_cond", "coverage_oracle", "mean_se", "se_ratio")
out <- out[order(out$sim, match(out$estimator, EST_ORDER), match(out$measure, MEAS), out$truth, out$t, out$tprime), ]
rownames(out) <- NULL
out$value <- signif(out$value, 8); out$mcse <- signif(out$mcse, 6); out$tau <- signif(out$tau, 8)
write.csv(out, "results/rebuild2_summary_long.csv", row.names = FALSE, na = "")
cat(nrow(out), "rows written to results/rebuild2_summary_long.csv\n")
print(table(out$measure, out$sim))
