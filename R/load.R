## Sources all code in R/. Usage, from the repository root: source("R/load.R")
local({
  dir <- if (dir.exists("R")) "R" else stop("run from the repository root")
  for (f in c("config.R", "dgp.R", "mlogit.R", "estimators.R", "variance.R", "truth.R", "rng.R", "run.R",
              "performance.R"))
    if (file.exists(file.path(dir, f))) sys.source(file.path(dir, f), envir = globalenv())
})
