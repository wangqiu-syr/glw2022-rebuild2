## Runs (or resumes) replications and writes per-replication files.
## Usage, from the repository root:
##   Rscript scripts/run_simulations.R <sims> <first_rep> <last_rep> <outdir> [cores]
## e.g. Rscript scripts/run_simulations.R sim1,sim2 1 20 results/pilot/runs 4
## Each replication is saved to <outdir>/<sim>/rep_NNNN.rds as it finishes; rerunning
## the same command skips finished replications. Progress goes to <outdir>/<sim>/run.log.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: run_simulations.R <sims> <first_rep> <last_rep> <outdir> [cores]")
source("R/load.R")
sims <- strsplit(args[1], ",")[[1]]
reps <- seq(as.integer(args[2]), as.integer(args[3]))
outdir <- args[4]
cores <- if (length(args) >= 5) as.integer(args[5]) else 4L
for (s in sims) {
  d <- file.path(outdir, s); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  run_simulation(SIMS[[s]], reps, d, cores = cores, log_file = file.path(d, "run.log"))
  cat(sprintf("%s: replications %d-%d finished in %.1f s wall-clock\n", s, min(reps), max(reps),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
