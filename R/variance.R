## Variance estimators for the pairwise effects (Stage 4; derivation in
## docs/decisions.md D26).
##
## A matching estimator is linear in the outcomes:
##   tau_hat(t,t') = (1/n) [ sum_{j in G_t} w_j Y_j - sum_{j in G_t'} w_j Y_j ],
## with w_j = number of units (including j itself) whose imputed Y(t) uses j
## (fractional under ties). Writing Y_j = m(V_j) + e_j, where V_j is what the
## estimator matches on (all covariates for CM, the scalar score for GPS), gives
##   tau_hat - tau = [ (1/n) sum_i tau(V_i) - tau ] + matching bias + (1/n) sum_j c_j e_j.
## The two random terms are uncorrelated, so
##   var_noise = (1/n^2) [ sum_{G_t} w_j^2 sig2_j + sum_{G_t'} w_j^2 sig2_j ],
##   var_sate  = (1/n^2) sum_s n_s Var(tau(V) | T = s)       (fixed group sizes, D3),
## where sig2_j is estimated by within-group nearest-neighbour differencing in
## the estimator's own metric, and Var(tau(V) | T = s) by the within-group
## variance of the unit-level imputed effects minus their average noise variance.
## var_noise alone is the variance given the covariates (for the sample-average
## effect); var_noise + var_sate is the variance for the population effect.
## Matching bias is not part of the variance; uncertainty from estimating the
## propensity score is ignored.

pair_variances_matching <- function(fit, T, A) {
  n <- length(T); pr <- t(combn(A, 2))
  V <- sapply(seq_len(A), function(t) fit$details[[t]]$v)        # n x A noise var of Yhat_i(t)
  wsig <- vapply(seq_len(A), function(t) {
    d <- fit$details[[t]]; sum(d$w^2 * d$sig2) }, 0)
  ns <- tabulate(T, A)
  out <- t(apply(pr, 1, function(q) {
    tau_i <- fit$Yhat[, q[1]] - fit$Yhat[, q[2]]
    v_i <- V[, q[1]] + V[, q[2]]
    var_noise <- (wsig[q[1]] + wsig[q[2]]) / n^2
    het <- vapply(seq_len(A), function(s) {
      g <- T == s; max(0, var(tau_i[g]) - mean(v_i[g])) }, 0)
    c(var_noise = var_noise, var_sate = sum(ns * het) / n^2)
  }))
  data.frame(t = pr[, 1], tprime = pr[, 2], var_noise = out[, 1], var_sate = out[, 2])
}

## DIF: independent groups of fixed size, so Var = s_t^2/n_t + s_t'^2/n_t'.
pair_variances_dif <- function(Y, T, A) {
  pr <- t(combn(A, 2))
  s2n <- vapply(seq_len(A), function(t) var(Y[T == t]) / sum(T == t), 0)
  data.frame(t = pr[, 1], tprime = pr[, 2], var_noise = NA_real_,
             var_sate = NA_real_, var_full = s2n[pr[, 1]] + s2n[pr[, 2]])
}
