## Estimators (docs/spec.md section 3). Each returns per-treatment means m_t;
## the pairwise estimate is tau_hat(t, t') = m_t - m_t'.

## DIF: sample difference (spec 3.1).
est_dif <- function(Y, T, A) vapply(seq_len(A), function(t) mean(Y[T == t]), 0)

## Nearest-neighbour imputation from a distance matrix D (n x m, rows = units to
## impute, columns = candidate donors with outcomes Yg). One neighbour; exact
## ties at the minimum (up to floating-point tolerance) are averaged (D10).
impute_nn <- function(D, Yg, rel_tol = 1e-10) {
  dmin <- apply(D, 1, min)
  tie <- D <= dmin + rel_tol * pmax(dmin, 1)    # recycles row-wise (column-major)
  as.vector((tie %*% Yg) / rowSums(tie))
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

## CM: covariates matching (spec 3.2). Distance = Mahalanobis on the X block
## + Mahalanobis on the Z block, covariances from the pooled sample (D9).
## Imputes Y_i(t) for every unit i and level t; averages over all n (D11).
est_cm <- function(W, T, Y, A, x_cols = c("x1", "x2", "x3", "x4"),
                   z_cols = c("z1", "z2", "z3")) {
  BX <- mahalanobis_whiten(W[, x_cols, drop = FALSE])
  BZ <- mahalanobis_whiten(W[, z_cols, drop = FALSE])
  Yhat <- matrix(NA_real_, length(Y), A)
  for (t in seq_len(A)) {
    idx <- which(T == t)
    D <- pairwise_dist(BX, idx) + pairwise_dist(BZ, idx)
    Yhat[, t] <- impute_nn(D, Y[idx])
  }
  list(m = colMeans(Yhat), Yhat = Yhat)
}

## GPS matching (spec 3.3, D13). S: n x A matrix of (estimated or true)
## propensity scores; for target t, match on the scalar S[, t].
est_gps <- function(S, T, Y, A) {
  Yhat <- matrix(NA_real_, length(Y), A)
  for (t in seq_len(A)) {
    idx <- which(T == t)
    D <- abs(outer(S[, t], S[idx, t], "-"))
    Yhat[, t] <- impute_nn(D, Y[idx])
  }
  list(m = colMeans(Yhat), Yhat = Yhat)
}

## Pairwise effects from per-treatment means.
pairs_from_m <- function(m) {
  A <- length(m); pr <- t(combn(A, 2))
  data.frame(t = pr[, 1], tprime = pr[, 2], est = m[pr[, 1]] - m[pr[, 2]])
}
