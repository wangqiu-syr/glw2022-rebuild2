## One replication, and a resumable parallel driver (spec 6, D17).

## Runs every estimator on one simulated data set. `seed` is the replication's
## L'Ecuyer-CMRG stream. Returns per-treatment means m (estimators x A) and
## per-replication diagnostics.
run_replication <- function(cfg, r, seed) {
  use_seed(seed)
  t0 <- proc.time()[["elapsed"]]
  smp <- simulate_sample(cfg)
  A <- cfg$A; W <- smp$W; T <- smp$T; Y <- smp$Y
  m <- matrix(NA_real_, length(ESTIMATORS), A, dimnames = list(ESTIMATORS, paste0("t", 1:A)))
  m["DIF", ] <- est_dif(Y, T, A)
  m["CM", ]  <- est_cm(W, T, Y, A)$m
  fit_main <- fit_mlogit(gps_design(W, "main"), T, A)
  m["GPS", ] <- est_gps(fit_main$fitted, T, Y, A)$m
  fit_cs <- fit_mlogit(gps_design(W, "cs"), T, A)
  m["GPS_cs", ] <- est_gps(fit_cs$fitted, T, Y, A)$m
  m["GPS_true", ] <- est_gps(smp$P_true, T, Y, A)$m
  list(sim = cfg$name, rep = r, m = m,
       mu_bar = colMeans(smp$mu),                                  # for the sample-average effect
       var_y = vapply(1:A, function(t) var(Y[T == t]), 0),          # for the D15 diagnostic
       diag = data.frame(sim = cfg$name, rep = r, n_candidates = smp$n_candidates,
                         main_converged = fit_main$converged, main_iter = fit_main$iter,
                         cs_converged = fit_cs$converged, cs_iter = fit_cs$iter,
                         min_phat_main = min(fit_main$fitted), min_phat_cs = min(fit_cs$fitted),
                         seconds = proc.time()[["elapsed"]] - t0))
}

rep_file <- function(outdir, r) file.path(outdir, sprintf("rep_%04d.rds", r))

## Runs replications `reps` of simulation `cfg`, in parallel, saving each
## replication to its own file as soon as it finishes. Replications whose
## file already exists are skipped, so an interrupted run can be resumed.
run_simulation <- function(cfg, reps, outdir, cores = 4L, log_file = NULL) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  seeds <- replication_seeds(cfg$master_seed, max(reps))
  todo <- reps[!file.exists(rep_file(outdir, reps))]
  logmsg <- function(...) {
    msg <- sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...))
    if (is.null(log_file)) cat(msg) else cat(msg, file = log_file, append = TRUE)
  }
  logmsg(cfg$name, ": ", length(reps) - length(todo), " done, ", length(todo), " to run")
  one <- function(r) {
    res <- run_replication(cfg, r, seeds[[r]])
    tmp <- paste0(rep_file(outdir, r), ".tmp")
    saveRDS(res, tmp); file.rename(tmp, rep_file(outdir, r))   # atomic write
    logmsg(cfg$name, " rep ", r, " (", round(res$diag$seconds, 2), " s)")
    r
  }
  out <- parallel::mclapply(todo, function(r) tryCatch(one(r), error = function(e) {
    logmsg(cfg$name, " rep ", r, " FAILED: ", conditionMessage(e)); NA_integer_ }),
    mc.cores = cores, mc.preschedule = FALSE)
  invisible(unlist(out))
}

## Reads saved replications into long data frames.
collect_simulation <- function(outdir) {
  files <- sort(list.files(outdir, pattern = "^rep_\\d+\\.rds$", full.names = TRUE))
  res <- lapply(files, readRDS)
  est <- do.call(rbind, lapply(res, function(x) {
    do.call(rbind, lapply(rownames(x$m), function(e)
      cbind(sim = x$sim, rep = x$rep, estimator = e, pairs_from_m(x$m[e, ]))))
  }))
  sate <- do.call(rbind, lapply(res, function(x) {
    p <- pairs_from_m(x$mu_bar); names(p)[3] <- "sate"; cbind(sim = x$sim, rep = x$rep, p)
  }))
  vary <- do.call(rbind, lapply(res, function(x) {
    A <- length(x$var_y); pr <- t(combn(A, 2))
    data.frame(sim = x$sim, rep = x$rep, t = pr[, 1], tprime = pr[, 2],
               varsum_sample = x$var_y[pr[, 1]] + x$var_y[pr[, 2]])
  }))
  diag <- do.call(rbind, lapply(res, `[[`, "diag"))
  list(estimates = est, sate = sate, varsum = vary, diag = diag)
}
