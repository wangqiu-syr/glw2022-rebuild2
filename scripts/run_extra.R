## Extra runs on the Stage-3 data sets (same random streams).
## Usage, from the repository root:
##   Rscript scripts/run_extra.R <mode> <sims> <first_rep> <last_rep> <outdir> [cores]
##   mode: cm_metrics (CM distance-metric sensitivity, D9) or stage4 (variances, D26)
## Resumable like run_simulations.R.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("usage: run_extra.R <mode> <sims> <first_rep> <last_rep> <outdir> [cores]")
source("R/load.R")
fun <- switch(args[1], cm_metrics = run_replication_cm_metrics, stage4 = run_replication_stage4,
              stop("unknown mode"))
reps <- seq(as.integer(args[3]), as.integer(args[4]))
cores <- if (length(args) >= 6) as.integer(args[6]) else 4L
for (s in strsplit(args[2], ",")[[1]]) {
  d <- file.path(args[5], s); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  run_simulation(SIMS[[s]], reps, d, cores = cores, log_file = file.path(d, "run.log"), rep_fun = fun)
  cat(sprintf("%s %s: replications %d-%d finished in %.1f s wall-clock\n", args[1], s, min(reps), max(reps),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
