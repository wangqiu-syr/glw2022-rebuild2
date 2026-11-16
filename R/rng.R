## Reproducible random-number streams, one per replication (spec 6, D17).
## Stream r = parallel::nextRNGStream applied r times to the L'Ecuyer-CMRG
## state obtained from set.seed(master_seed).
replication_seeds <- function(master_seed, n) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(master_seed)
  s <- .Random.seed; out <- vector("list", n)
  for (r in seq_len(n)) { s <- parallel::nextRNGStream(s); out[[r]] <- s }
  out
}

use_seed <- function(seed) {
  RNGkind("L'Ecuyer-CMRG")
  assign(".Random.seed", seed, envir = globalenv())
}
