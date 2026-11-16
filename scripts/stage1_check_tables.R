## Stage 1: internal-consistency check of the transcribed tables.
## For DIF, CM and GPS, tau_hat(t,t') = m_t - m_t' where m_t depends only on t
## (the imputed Y_i(t) does not depend on t'), so E(tau_hat(t,t')) - tau(t,t')
## must be exactly of the form h_t - h_t'.  A transcription error in a single
## bias cell would show up as a large residual.  Run from the repository root.
t3 <- read.csv("data/paper/table3_sim2.csv"); t1 <- read.csv("data/paper/table1_sim1.csv")
additive_fit <- function(t, tp, y, A) {
  X <- matrix(0, length(y), A); X[cbind(seq_along(y), t)] <- 1; X[cbind(seq_along(y), tp)] <- -1
  f <- lm.fit(X[, -1, drop = FALSE], y); r <- f$residuals
  list(h = c(0, f$coefficients), max_abs_resid = max(abs(r)), worst = which.max(abs(r)))
}
cat("Table 3 (Sim 2): additivity of signed bias columns, 36 pairs, 8 free parameters\n")
for (col in c("dif_bias", "cm_bias", "gps_bias")) {
  f <- additive_fit(t3$t, t3$tprime, t3[[col]], 9)
  cat(sprintf("  %-8s max |residual| = %.5f (worst: %s)\n", col, f$max_abs_resid, t3$effect[f$worst]))
}
## Table 1 reports |bias| only.  Search all 2^6 sign patterns for an additive fit.
cat("Table 1 (Sim 1): best sign pattern making |bias| additive (6 pairs, 3 free parameters)\n")
S <- as.matrix(expand.grid(rep(list(c(-1, 1)), 6)))
for (col in c("dif_absbias", "cm_absbias", "gps_absbias")) {
  res <- apply(S, 1, function(s) additive_fit(t1$t, t1$tprime, s * t1[[col]], 4)$max_abs_resid)
  best <- which.min(res)
  cat(sprintf("  %-11s min over sign patterns of max |residual| = %.5f; signs (up to global flip): %s\n",
              col, res[best], paste(S[best, ] * S[best, 1], collapse = " ")))
}

## Save the per-treatment decomposition h_t (h_1 = 0) of the paper's Table 3 biases.
dec <- data.frame(t = 1:9, sapply(c(dif = "dif_bias", cm = "cm_bias", gps = "gps_bias"),
                  function(col) round(additive_fit(t3$t, t3$tprime, t3[[col]], 9)$h, 4)))
names(dec)[-1] <- paste0("h_", names(dec)[-1], "_rel_t1")
write.csv(dec, "results/stage1/paper_table3_bias_decomposition.csv", row.names = FALSE)
