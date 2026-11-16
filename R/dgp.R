## Data-generating process (docs/spec.md section 2).

COV_NAMES <- c("x1", "x2", "x3", "x4", "z1", "z2", "z3")

SIGMA_X <- matrix(c( 2.0,  1.0, -1.0,  0.5,
                     1.0,  1.0, -0.5,  0.2,
                    -1.0, -0.5,  1.0, -0.3,
                     0.5,  0.2, -0.3,  2.0), 4, 4, byrow = TRUE)
CHOL_SIGMA_X <- chol(SIGMA_X)          # upper triangular R with R'R = Sigma_X

## E[W] under spec 2.1: X mean 0, Z1 ~ U(-3,3), Z2 ~ chi2_1, Z3 ~ Bern(0.5)
E_W <- c(x1 = 0, x2 = 0, x3 = 0, x4 = 0, z1 = 0, z2 = 1, z3 = 0.5)

## Squared terms in eta_t (paper p. 974). eta_1 = 0; for t >= 2,
## eta_t = 0.1 * sum over all seven covariates, with the listed ones squared.
ETA_SQUARED <- list(
  sim1 = list(NULL, "x1", c("x2", "x3"), c("x2", "x3", "x4")),
  sim2 = list(NULL, "x1", c("x2", "x3"), c("x2", "x3", "x4"),
              c("x2", "x3", "x4", "z1"),
              c("x2", "x3", "x4", "z1", "z2"),
              c("x2", "x3", "x4", "z1", "z2", "z3"),
              COV_NAMES,
              character(0))
)

## Outcome models mu_t(w) = intercept_t + coef_t' w (paper pp. 974-975).
## Sim 2 is written in z2 - 1 and z3 - 0.5, i.e. intercept -(b_z2 + 0.5 b_z3).
OUTCOME_COEF <- list(
  sim1 = rbind(c(1, 1, 1, 1, 1, 1, 1),
               c(2, 3, 1, 2, 2, 2, 2),
               c(3, 1, 2, 1, 3, 1, 3),
               c(4, 1, 3, 4, 4, 1, 4)),
  sim2 = rbind(c(1, 1, 1, 1,   1, 1, 1),
               c(2, 3, 1, 2,   2, 2, 2),
               c(3, 1, 2, 1,   3, 1, 3),
               c(4, 1, 3, 4,   4, 1, 4),
               c(5, 1, 1, 1,   1, 1, 1),
               c(6, 3, 1, 2,   2, 2, 9),
               c(7, 1, 8, 1,   3, 1, 8),
               c(8, 1, 5, 2.5, 4, 1, 7),
               c(9, 1, 7, 1.5, 4, 1, 7))
)
OUTCOME_COEF <- lapply(OUTCOME_COEF, function(B) { colnames(B) <- COV_NAMES; B })
OUTCOME_INTERCEPT <- list(
  sim1 = rep(0, 4),
  sim2 = -(OUTCOME_COEF$sim2[, "z2"] * 1 + OUTCOME_COEF$sim2[, "z3"] * 0.5)
)

## Draw N i.i.d. covariate vectors W = (X, Z). Draw order: X (4N normals), Z1, Z2, Z3.
draw_covariates <- function(N) {
  X  <- matrix(rnorm(4 * N), N, 4) %*% CHOL_SIGMA_X
  z1 <- runif(N, -3, 3)
  z2 <- rchisq(N, df = 1)
  z3 <- as.numeric(rbinom(N, 1, 0.5))
  W <- cbind(X, z1, z2, z3)
  colnames(W) <- COV_NAMES
  W
}

## N x A matrix of eta_t(W).
eta_matrix <- function(W, sim) {
  sq <- ETA_SQUARED[[sim]]
  A <- length(sq)
  E <- matrix(0, nrow(W), A)
  for (t in seq_len(A)[-1]) {
    V <- W[, COV_NAMES, drop = FALSE]
    if (length(sq[[t]])) V[, sq[[t]]] <- V[, sq[[t]]]^2
    E[, t] <- 0.1 * rowSums(V)
  }
  E
}

## Softmax of a matrix of linear predictors, row-wise, numerically stable.
softmax_rows <- function(E) {
  E <- E - E[cbind(seq_len(nrow(E)), max.col(E, ties.method = "first"))]
  P <- exp(E)
  P / rowSums(P)
}

## True generalized propensity scores p_t(W), N x A (D2: standard softmax, eta_1 = 0).
true_gps <- function(W, sim) softmax_rows(eta_matrix(W, sim))

## N x A matrix of mu_t(W) = E[Y(t) | W].
mu_matrix <- function(W, sim) {
  M <- W[, COV_NAMES, drop = FALSE] %*% t(OUTCOME_COEF[[sim]])
  sweep(M, 2, OUTCOME_INTERCEPT[[sim]], "+")
}

## Categorical draw by inverse CDF: T_i = min{t : u_i <= P_i1 + ... + P_it}.
draw_treatment <- function(P, u) {
  A <- ncol(P)
  if (A == 1L) return(rep(1L, nrow(P)))
  cumP <- P[, -A, drop = FALSE]
  if (A > 2L) for (k in 2:(A - 1L)) cumP[, k] <- cumP[, k - 1L] + P[, k]
  as.integer(rowSums(u > cumP)) + 1L
}

## One simulated data set (spec 2.4, D3): candidates (W, T) are drawn in batches
## of `batch_size` from the covariate distribution and the assignment model; a
## candidate is kept if its group still has fewer than n_t members. Stops when
## every group is full. Then eps ~ N(0,1) is drawn for the kept units, in order.
## Returns units sorted by treatment, then by order of drawing.
simulate_sample <- function(cfg, batch_size = BATCH_SIZE) {
  sim <- cfg$name; A <- cfg$A; nt <- cfg$n_per_group
  kept_W <- vector("list", A); count <- integer(A); n_cand <- 0L
  while (any(count < nt)) {
    W <- draw_covariates(batch_size)
    u <- runif(batch_size)
    T <- draw_treatment(true_gps(W, sim), u)
    n_cand <- n_cand + batch_size
    for (t in which(count < nt)) {
      idx <- which(T == t)
      take <- idx[seq_len(min(length(idx), nt - count[t]))]
      if (length(take)) {
        kept_W[[t]] <- rbind(kept_W[[t]], W[take, , drop = FALSE])
        count[t] <- count[t] + length(take)
      }
    }
  }
  W <- do.call(rbind, kept_W)
  T <- rep(seq_len(A), each = nt)
  mu <- mu_matrix(W, sim)
  eps <- rnorm(nrow(W))
  Y <- mu[cbind(seq_len(nrow(W)), T)] + eps
  list(W = W, T = T, Y = Y, mu = mu, P_true = true_gps(W, sim), n_candidates = n_cand)
}
