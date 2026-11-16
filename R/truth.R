## True effects (docs/spec.md section 4).

## tau_super: population effect E[Y(t) - Y(t')] under spec 2.1, analytic
## (mu_t is linear in W and E[W] is known).
tau_super <- function(sim) {
  EY <- OUTCOME_INTERCEPT[[sim]] + as.vector(OUTCOME_COEF[[sim]] %*% E_W)
  p <- pairs_from_m(EY); names(p)[3] <- "tau_super"; p
}

## Population moments by Monte Carlo integration, weighting by the true p_t
## instead of sampling T: E[g(W) | T = t] = E[g(W) p_t(W)] / E[p_t(W)].
## Uses one L'Ecuyer-CMRG stream per chunk; chunks run in parallel.
## Returns pi_t, E[W | T = t] (A x 7) and E[W W' | T = t] (A x 7 x 7) with MC SEs
## of the conditional means.
population_moments <- function(sim, n_chunks = 100L, chunk_size = 1e6, seed = 1L,
                               cores = 4L) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  seeds <- vector("list", n_chunks); s <- .Random.seed
  for (c in seq_len(n_chunks)) { s <- parallel::nextRNGStream(s); seeds[[c]] <- s }
  one <- function(c) {
    assign(".Random.seed", seeds[[c]], envir = globalenv())
    W <- draw_covariates(chunk_size); P <- true_gps(W, sim)
    Wc <- sweep(W, 2, E_W)                     # centring reduces MC variance
    list(sp = colSums(P), swp = crossprod(P, Wc),
         swwp = lapply(seq_len(ncol(P)), function(t) crossprod(Wc * P[, t], Wc)),
         ## per-chunk conditional means, for a between-chunk standard error
         cm = sweep(crossprod(P, Wc), 1, colSums(P), "/"))
  }
  res <- parallel::mclapply(seq_len(n_chunks), one, mc.cores = cores)
  A <- length(ETA_SQUARED[[sim]]); N <- n_chunks * chunk_size
  sp <- Reduce(`+`, lapply(res, `[[`, "sp")); swp <- Reduce(`+`, lapply(res, `[[`, "swp"))
  EWc_T <- sweep(swp, 1, sp, "/")              # E[W - E W | T = t]
  EW_T <- sweep(EWc_T, 2, E_W, "+")
  cm <- simplify2array(lapply(res, `[[`, "cm"))  # A x 7 x chunks
  se_EW_T <- apply(cm, c(1, 2), sd) / sqrt(n_chunks)
  covW_T <- lapply(seq_len(A), function(t) {
    S <- Reduce(`+`, lapply(res, function(r) r$swwp[[t]])) / sp[t]
    S - tcrossprod(EWc_T[t, ])                 # Cov(W | T = t)
  })
  dimnames(EW_T) <- dimnames(se_EW_T) <- list(paste0("t", seq_len(A)), COV_NAMES)
  cm_chunks <- sweep(cm, 2, E_W, "+")          # per-chunk E[W | T = t]
  list(sim = sim, N = N, pi = sp / N, EW_T = EW_T, se_EW_T = se_EW_T, covW_T = covW_T,
       EW_T_chunks = cm_chunks)
}

## Derived population quantities for all pairs:
##   tau_super, tau_mix (equal-weight mixture of the groups, spec 4),
##   E_dif = E[Ybar_t - Ybar_t'] (population DIF expectation),
##   varY_T = Var(Y | T = t) incl. eps (for the D15 diagnostic),
##   MC SEs of tau_mix and E_dif.
population_truth <- function(mom) {
  sim <- mom$sim; Bc <- OUTCOME_COEF[[sim]]; a <- OUTCOME_INTERCEPT[[sim]]
  A <- nrow(Bc)
  m_mix <- colMeans(mom$EW_T)
  EY_mix <- a + as.vector(Bc %*% m_mix)
  EY_own <- a + rowSums(Bc * mom$EW_T)         # E[mu_t(W) | T = t]
  varY_T <- vapply(seq_len(A), function(t) as.numeric(Bc[t, ] %*% mom$covW_T[[t]] %*% Bc[t, ]) + 1, 0)
  ## Monte Carlo SEs from the spread of the same quantities across chunks
  nc <- dim(mom$EW_T_chunks)[3]
  per_chunk <- lapply(seq_len(nc), function(c) {
    E <- mom$EW_T_chunks[, , c]
    list(mix = a + as.vector(Bc %*% colMeans(E)), own = a + rowSums(Bc * E))
  })
  pr <- t(combn(A, 2))
  chunk_se <- function(f) apply(sapply(per_chunk, function(pc) f(pc)), 1, sd) / sqrt(nc)
  out <- tau_super(sim)
  out$tau_mix <- EY_mix[pr[, 1]] - EY_mix[pr[, 2]]
  out$se_tau_mix <- chunk_se(function(pc) pc$mix[pr[, 1]] - pc$mix[pr[, 2]])
  out$E_dif <- EY_own[pr[, 1]] - EY_own[pr[, 2]]
  out$se_E_dif <- chunk_se(function(pc) pc$own[pr[, 1]] - pc$own[pr[, 2]])
  out$dif_bias_pop <- out$E_dif - out$tau_super
  out$varY_t <- varY_T[pr[, 1]]; out$varY_tprime <- varY_T[pr[, 2]]
  out$varsum <- out$varY_t + out$varY_tprime
  cbind(sim = sim, out)
}
