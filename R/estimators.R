## Estimators (docs/spec.md section 3). Each returns per-treatment means m_t;
## the pairwise estimate is tau_hat(t, t') = m_t - m_t'.

## DIF: sample difference (spec 3.1).
est_dif <- function(Y, T, A) vapply(seq_len(A), function(t) mean(Y[T == t]), 0)

## Nearest-neighbour imputation from a distance matrix D (n x m, rows = units to
## impute, columns = candidate donors with outcomes Yg). One neighbour; exact
## ties at the minimum (up to floating-point tolerance) are averaged (D10).
impute_nn <- function(D, Yg, rel_tol = 1e-10) as.vector(match_matrix(D, rel_tol) %*% Yg)

## Row-normalised match matrix: a[i, j] = 1/k if donor j is one of the k donors
## tied at unit i's minimum distance, else 0.
match_matrix <- function(D, rel_tol = 1e-10) {
  dmin <- apply(D, 1, min)
  tie <- D <= dmin + rel_tol * pmax(dmin, 1)    # recycles row-wise (column-major)
  tie / rowSums(tie)
}

## Matching details used by the variance estimator (R/variance.R, Stage 4).
## D: n x n_t distances from all units to the donors idx (group t); Y: all outcomes.
## Returns donor weights w_j = sum_i a_ij, and within-group nearest-neighbour
## estimates sig2_j = (Y_j - Y_l)^2 / 2 (l = j's nearest other donor, same metric).
match_details <- function(D, idx, Y, rel_tol = 1e-10) {
  a <- match_matrix(D, rel_tol)
  Dw <- D[idx, , drop = FALSE]; diag(Dw) <- Inf
  aw <- match_matrix(Dw, rel_tol)
  Yg <- Y[idx]
  sig2 <- rowSums(aw * outer(Yg, Yg, "-")^2) / 2
  list(idx = idx, w = colSums(a), sig2 = sig2, v = as.vector(a^2 %*% sig2))
}

## Rows of M multiplied so that Euclidean distance = Mahalanobis distance with
## V = (1/n) sum (M_i - Mbar)(M_i - Mbar)' over all rows (spec 3.2, D9).
mahalanobis_whiten <- function(M) {
  M <- as.matrix(M)
  V <- crossprod(sweep(M, 2, colMeans(M))) / nrow(M)
  M %*% solve(chol(V))
}

## Euclidean distances between rows of B (n x d) and rows B[idx, ] (n x |idx|),
## computed from coordinate differences, so a unit's distance to itself is exactly 0.
pairwise_dist <- function(B, idx) {
  D2 <- 0
  for (k in seq_len(ncol(B))) D2 <- D2 + outer(B[, k], B[idx, k], "-")^2
  sqrt(D2)
}

## CM: covariates matching (spec 3.2). Covariances from the pooled sample (D9).
## metric = "blocks" (primary, D9a): Mahalanobis on X + Mahalanobis on Z;
##          "stacked" (D9b): one Mahalanobis distance on (X, Z);
##          "x_mahal_z_euclid" (D9c): Mahalanobis on X + Euclidean on Z.
## Imputes Y_i(t) for every unit i and level t; averages over all n (D11).
## details = TRUE adds what the Stage-4 variance estimator needs.
cm_blocks <- function(W, metric, x_cols, z_cols) {
  switch(metric,
         blocks = list(mahalanobis_whiten(W[, x_cols, drop = FALSE]),
                       mahalanobis_whiten(W[, z_cols, drop = FALSE])),
         stacked = list(mahalanobis_whiten(W[, c(x_cols, z_cols), drop = FALSE])),
         x_mahal_z_euclid = list(mahalanobis_whiten(W[, x_cols, drop = FALSE]),
                                 as.matrix(W[, z_cols, drop = FALSE])),
         stop("unknown metric"))
}
cm_dist <- function(blocks, idx) Reduce(`+`, lapply(blocks, pairwise_dist, idx = idx))

est_cm <- function(W, T, Y, A, x_cols = c("x1", "x2", "x3", "x4"),
                   z_cols = c("z1", "z2", "z3"), metric = "blocks", details = FALSE) {
  B <- cm_blocks(W, metric, x_cols, z_cols)
  Yhat <- matrix(NA_real_, length(Y), A); det <- vector("list", A)
  for (t in seq_len(A)) {
    idx <- which(T == t)
    D <- cm_dist(B, idx)
    Yhat[, t] <- impute_nn(D, Y[idx])
    if (details) det[[t]] <- match_details(D, idx, Y)
  }
  out <- list(m = colMeans(Yhat), Yhat = Yhat)
  if (details) out$details <- det
  out
}

## GPS matching (spec 3.3, D13). S: n x A matrix of (estimated or true)
## propensity scores; for target t, match on the scalar S[, t].
est_gps <- function(S, T, Y, A, details = FALSE) {
  Yhat <- matrix(NA_real_, length(Y), A); det <- vector("list", A)
  for (t in seq_len(A)) {
    idx <- which(T == t)
    D <- abs(outer(S[, t], S[idx, t], "-"))
    Yhat[, t] <- impute_nn(D, Y[idx])
    if (details) det[[t]] <- match_details(D, idx, Y)
  }
  out <- list(m = colMeans(Yhat), Yhat = Yhat)
  if (details) out$details <- det
  out
}

## Pairwise effects from per-treatment means.
pairs_from_m <- function(m) {
  A <- length(m); pr <- t(combn(A, 2))
  data.frame(t = pr[, 1], tprime = pr[, 2], est = m[pr[, 1]] - m[pr[, 2]])
}
