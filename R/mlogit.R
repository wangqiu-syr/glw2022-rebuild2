## Multinomial logistic regression by maximum likelihood (Newton-Raphson).
## Used to estimate the generalized propensity score (spec 3.3, D12).
## Baseline category 1: log{p_k(x) / p_1(x)} = x' b_k, k = 2..A.

## Design matrices (without intercept; fit_mlogit adds it).
##   "main": main effects of the 7 covariates (primary GPS model, D12a)
##   "cs":   main effects + x1^2..x4^2, z1^2, z2^2 (correctly specified, D12b;
##           z3^2 = z3 is omitted because it is collinear)
gps_design <- function(W, type = c("main", "cs")) {
  type <- match.arg(type)
  D <- W[, COV_NAMES, drop = FALSE]
  if (type == "cs") {
    S <- W[, c("x1", "x2", "x3", "x4", "z1", "z2"), drop = FALSE]^2
    colnames(S) <- paste0(colnames(S), "_sq")
    D <- cbind(D, S)
  }
  D
}

## Fit. X: n x q design (no intercept column); T: integer in 1..A.
## Columns are standardised internally (fitted probabilities are unaffected).
## Returns fitted n x A probabilities, the log-likelihood, iterations and a
## convergence flag (max |Newton step| < tol, on the standardised scale).
fit_mlogit <- function(X, T, A, maxit = 100L, tol = 1e-8) {
  X <- as.matrix(X)
  ctr <- colMeans(X); scl <- apply(X, 2, sd); scl[scl == 0] <- 1
  Xs <- cbind(1, sweep(sweep(X, 2, ctr), 2, scl, "/"))
  n <- nrow(Xs); p <- ncol(Xs); K <- A - 1L
  Yind <- matrix(0, n, A); Yind[cbind(seq_len(n), T)] <- 1
  B <- matrix(0, p, K)
  probs <- function(B) softmax_rows(cbind(0, Xs %*% B))
  loglik <- function(P) sum(log(pmax(P[cbind(seq_len(n), T)], .Machine$double.xmin)))
  P <- probs(B); ll <- loglik(P); converged <- FALSE
  for (iter in seq_len(maxit)) {
    g <- as.vector(crossprod(Xs, Yind[, -1, drop = FALSE] - P[, -1, drop = FALSE]))
    H <- matrix(0, p * K, p * K)            # negative Hessian (information)
    for (k in seq_len(K)) for (l in k:K) {
      w <- P[, k + 1] * ((k == l) - P[, l + 1])
      blk <- crossprod(Xs, Xs * w)
      ik <- (k - 1) * p + seq_len(p); il <- (l - 1) * p + seq_len(p)
      H[ik, il] <- blk; H[il, ik] <- t(blk)
    }
    step <- tryCatch(solve(H, g), error = function(e) qr.solve(H, g, tol = 1e-12))
    ## step halving if the log-likelihood does not increase
    s <- 1
    repeat {
      B_new <- B + s * matrix(step, p, K)
      P_new <- probs(B_new); ll_new <- loglik(P_new)
      if (ll_new >= ll - 1e-10 || s < 1e-6) break
      s <- s / 2
    }
    B <- B_new; P <- P_new; ll <- ll_new
    if (max(abs(s * step)) < tol) { converged <- TRUE; break }
  }
  list(fitted = P, loglik = ll, iter = iter, converged = converged,
       coef_std = B, center = ctr, scale = scl)
}
