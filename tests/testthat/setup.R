## Load the project code. testthat runs tests with the working directory set to
## tests/testthat, so source relative to the repository root.
local({
  root <- normalizePath(file.path("..", ".."))
  old <- setwd(root); on.exit(setwd(old))
  source(file.path(root, "R", "load.R"))
})
