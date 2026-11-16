## Summarises the pilot run (timing, convergence, preview of bias).
## From the repository root, after run_simulations.R wrote results/pilot/runs:
##   Rscript scripts/pilot_summary.R
source("R/load.R")
out <- "results/pilot"
lines <- c("# Pilot run (Stage 2)", "",
           sprintf("R %s; %d cores; generated %s.", getRversion(), parallel::detectCores(), format(Sys.Date())), "")
timing <- list(); perf_all <- list()
for (s in names(SIMS)) {
  col <- collect_simulation(file.path(out, "runs", s))
  dg <- col$diag; R <- nrow(dg)
  tr <- read.csv(sprintf("results/truth/%s_truth.csv", s))
  p_super <- performance(col$estimates, tr, "tau_super")
  p_mix <- performance(col$estimates, tr, "tau_mix")
  perf <- merge(p_super, p_mix[, c("estimator", "t", "tprime", "bias")], by = c("estimator", "t", "tprime"),
                suffixes = c("", "_vs_mix"))
  perf <- merge(perf, tr[, c("t", "tprime", "dif_bias_pop")], by = c("t", "tprime"))
  perf_all[[s]] <- perf
  write.csv(perf, file.path(out, sprintf("pilot_performance_%s.csv", s)), row.names = FALSE)
  timing[[s]] <- data.frame(sim = s, reps = R, mean_sec_per_rep = mean(dg$seconds),
                            max_sec_per_rep = max(dg$seconds),
                            projected_500_reps_4cores_min = 500 * mean(dg$seconds) / 4 / 60,
                            mean_candidates = mean(dg$n_candidates),
                            main_logit_converged = sum(dg$main_converged),
                            cs_logit_converged = sum(dg$cs_converged),
                            max_iter_main = max(dg$main_iter), max_iter_cs = max(dg$cs_iter),
                            min_phat_cs = min(dg$min_phat_cs))
  ## DIF check against the population value
  d <- perf[perf$estimator == "DIF", ]
  z <- (d$bias - d$dif_bias_pop) / d$mcse_bias
  agg <- aggregate(cbind(abs_bias, mse, mcse_bias) ~ estimator, perf, mean)
  agg <- agg[match(ESTIMATORS, agg$estimator), ]
  lines <- c(lines, sprintf("## %s (%d replications)", s, R), "",
             "Averages over all pairs (bias against tau_super; with only 20 replications these are rough):", "",
             "| estimator | mean abs bias | mean MSE | mean MC SE of bias |", "|---|---|---|---|",
             sprintf("| %s | %.3f | %.3f | %.3f |", agg$estimator, agg$abs_bias, agg$mse, agg$mcse_bias), "",
             sprintf("DIF sanity check: simulated DIF bias minus population DIF bias, in MC-SE units: mean %.2f, max |z| %.2f over %d pairs.",
                     mean(z), max(abs(z)), length(z)), "")
}
tm <- do.call(rbind, timing)
write.csv(tm, file.path(out, "pilot_timing.csv"), row.names = FALSE)
lines <- c(lines, "## Timing and convergence", "",
           "| sim | reps | s/rep (mean) | s/rep (max) | projected 500 reps on 4 cores (min) | candidates drawn (mean) | main logit converged | cs logit converged | max iter main / cs | min p-hat (cs) |",
           "|---|---|---|---|---|---|---|---|---|---|",
           sprintf("| %s | %d | %.2f | %.2f | %.1f | %.0f | %d/%d | %d/%d | %d / %d | %.1e |", tm$sim, tm$reps,
                   tm$mean_sec_per_rep, tm$max_sec_per_rep, tm$projected_500_reps_4cores_min, tm$mean_candidates,
                   tm$main_logit_converged, tm$reps, tm$cs_logit_converged, tm$reps, tm$max_iter_main,
                   tm$max_iter_cs, tm$min_phat_cs))
writeLines(lines, file.path(out, "pilot_summary.md"))
cat(lines, sep = "\n")
