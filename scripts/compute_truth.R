## Population truths and reference quantities (docs/spec.md section 4):
## tau_super (analytic), tau_mix, population DIF expectation and Var(Y | T = t),
## by Monte Carlo integration with 1e8 covariate draws per simulation.
## From the repository root:  Rscript scripts/compute_truth.R
source("R/load.R")
dir.create("results/truth", recursive = TRUE, showWarnings = FALSE)
for (s in names(SIMS)) {
  t0 <- Sys.time()
  mom <- population_moments(s, n_chunks = 100L, chunk_size = 1e6, seed = 777L, cores = 4L)
  tr <- population_truth(mom)
  write.csv(tr, sprintf("results/truth/%s_truth.csv", s), row.names = FALSE)
  pt <- data.frame(t = seq_along(mom$pi), pi = mom$pi, round(mom$EW_T, 6))
  write.csv(pt, sprintf("results/truth/%s_group_moments.csv", s), row.names = FALSE)
  cat(s, "done in", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s; max SE(tau_mix) =",
      signif(max(tr$se_tau_mix), 2), "; max SE(E_dif) =", signif(max(tr$se_E_dif), 2), "\n")
}
