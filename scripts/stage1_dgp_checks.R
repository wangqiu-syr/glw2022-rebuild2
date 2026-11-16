## Stage 1 diagnostic: which reading of the paper's DGP reproduces the paper's
## own DIF (sample-difference) column?  The DIF estimator involves no modelling
## choices, so E[DIF] - tau depends only on the DGP and on the definition of the
## true effect tau.  We compute population quantities by Monte Carlo
## integration (large N, weighting by the true treatment probabilities rather
## than sampling T, to reduce noise).
##
## Run:  Rscript scripts/stage1_dgp_checks.R
suppressPackageStartupMessages(library(MASS))

Sigma_X <- matrix(c( 2.0,  1.0, -1.0,  0.5,
                     1.0,  1.0, -0.5,  0.2,
                    -1.0, -0.5,  1.0, -0.3,
                     0.5,  0.2, -0.3,  2.0), 4, 4, byrow = TRUE)

draw_W <- function(N, z1_type) {
  X <- MASS::mvrnorm(N, rep(0, 4), Sigma_X)
  z1 <- switch(z1_type,
               unif_m3_3   = runif(N, -3, 3),
               norm_m3_var3 = rnorm(N, -3, sqrt(3)),
               norm_m3_sd3  = rnorm(N, -3, 3),
               norm_0_var3  = rnorm(N, 0, sqrt(3)))
  z2 <- rchisq(N, 1)
  z3 <- rbinom(N, 1, 0.5)
  cbind(x1 = X[, 1], x2 = X[, 2], x3 = X[, 3], x4 = X[, 4], z1 = z1, z2 = z2, z3 = z3)
}

## linear predictors eta_t (paper p. 974); 'sq' marks squared terms
eta_spec <- list(
  sim1 = list(c(), c("x1"), c("x2", "x3"), c("x2", "x3", "x4")),
  sim2 = list(NULL, c("x1"), c("x2", "x3"), c("x2", "x3", "x4"),
              c("x2", "x3", "x4", "z1"), c("x2", "x3", "x4", "z1", "z2"),
              c("x2", "x3", "x4", "z1", "z2", "z3"),
              c("x1", "x2", "x3", "x4", "z1", "z2", "z3"), c()))
eta_mat <- function(W, sim) {
  sp <- eta_spec[[sim]]
  E <- matrix(0, nrow(W), length(sp))
  for (t in seq_along(sp)[-1]) {
    V <- W
    if (length(sp[[t]])) V[, sp[[t]]] <- V[, sp[[t]]]^2
    E[, t] <- 0.1 * rowSums(V)
  }
  E
}
probs <- function(W, sim) {          # softmax, eta_1 = 0 (baseline)
  E <- eta_mat(W, sim); E <- E - apply(E, 1, max); P <- exp(E); P / rowSums(P)
}

## outcome coefficient vectors on (x1..x4, z1, z2, z3) and intercepts
beta <- list(
  sim1 = rbind(c(1,1,1,1,1,1,1), c(2,3,1,2,2,2,2), c(3,1,2,1,3,1,3), c(4,1,3,4,4,1,4)),
  sim2 = rbind(c(1,1,1,1,1,1,1), c(2,3,1,2,2,2,2), c(3,1,2,1,3,1,3), c(4,1,3,4,4,1,4),
               c(5,1,1,1,1,1,1), c(6,3,1,2,2,2,9), c(7,1,8,1,3,1,8), c(8,1,5,2.5,4,1,7),
               c(9,1,7,1.5,4,1,7)))
## Sim 2 writes z2 - 1 and z3 - 0.5, i.e. intercept -(b_z2 * 1 + b_z3 * 0.5)
intercept <- function(sim) if (sim == "sim1") rep(0, 4) else -(beta$sim2[, 6] + 0.5 * beta$sim2[, 7])
mu <- function(W, sim) sweep(W %*% t(beta[[sim]]), 2, intercept(sim), "+")

pop_summary <- function(sim, z1_type, N = 4e6, seed = 1) {
  set.seed(seed)
  W <- draw_W(N, z1_type); P <- probs(W, sim); M <- mu(W, sim)
  A <- ncol(P); pi_t <- colMeans(P)
  ## E[mu_s(W) | T = t] for all s, t  (A x A: rows t = group, cols s = outcome)
  cond_mu  <- t(sapply(1:A, function(t) colSums(M * P[, t]) / sum(P[, t])))
  cond_mu2 <- t(sapply(1:A, function(t) colSums(M^2 * P[, t]) / sum(P[, t])))
  cond_var_Y <- diag(cond_mu2) - diag(cond_mu)^2 + 1        # Var(Y | T = t), incl. eps
  super_mu <- colMeans(M)                                    # E[mu_s(W)]
  mix_mu   <- colMeans(cond_mu)                              # equal-weight mixture of groups
  pr <- t(combn(A, 2))
  data.frame(sim = sim, z1 = z1_type, t = pr[, 1], tprime = pr[, 2],
             tau_super = super_mu[pr[, 1]] - super_mu[pr[, 2]],
             tau_mix   = mix_mu[pr[, 1]] - mix_mu[pr[, 2]],
             E_dif     = diag(cond_mu)[pr[, 1]] - diag(cond_mu)[pr[, 2]],
             varsum    = cond_var_Y[pr[, 1]] + cond_var_Y[pr[, 2]],
             pi_t = paste(round(pi_t, 3), collapse = "/"))
}

z1_types <- c("unif_m3_3", "norm_m3_var3", "norm_m3_sd3", "norm_0_var3")
res <- do.call(rbind, lapply(c("sim1", "sim2"), function(s)
  do.call(rbind, lapply(z1_types, function(z) pop_summary(s, z)))))
res$bias_super <- res$E_dif - res$tau_super
res$bias_mix   <- res$E_dif - res$tau_mix
write.csv(res, "results/stage1/dgp_population_checks.csv", row.names = FALSE)

t1 <- read.csv("data/paper/table1_sim1.csv"); t3 <- read.csv("data/paper/table3_sim2.csv")
cmp <- function(sim, z, truth) {
  r <- res[res$sim == sim & res$z1 == z, ]
  b <- r[[paste0("bias_", truth)]]
  if (sim == "sim1") { p <- t1$dif_absbias; b <- abs(b) } else p <- t3$dif_bias
  data.frame(sim = sim, z1 = z, truth = truth, n_effects = length(p),
             rmsd = sqrt(mean((b - p)^2)), max_absdiff = max(abs(b - p)),
             cor = if (sd(b) > 0) cor(b, p) else NA,
             rmsd_signflip = if (sim == "sim2") sqrt(mean((-b - p)^2)) else NA,
             pi_t = r$pi_t[1])
}
summ <- do.call(rbind, lapply(c("sim1", "sim2"), function(s)
  do.call(rbind, lapply(z1_types, function(z)
    do.call(rbind, lapply(c("super", "mix"), function(tr) cmp(s, z, tr)))))))
write.csv(summ, "results/stage1/dgp_population_checks_summary.csv", row.names = FALSE)
print(summ, digits = 3, row.names = FALSE)
