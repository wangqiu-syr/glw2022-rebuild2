## Stage 1: evidence on the sampling design and the true effect.
## If each treatment group is sampled to exactly n_t units (design B) the pooled
## sample is an equal-weight mixture of the groups, so matching estimators that
## average over the whole sample target tau_mix; if the paper measured bias
## against the population effect tau_super, its CM biases should carry the
## signature tau_mix - tau_super.  Uses results/stage1/dgp_population_checks.csv
## (from stage1_dgp_checks.R).  Run from the repository root.
r  <- read.csv("results/stage1/dgp_population_checks.csv")
t3 <- read.csv("data/paper/table3_sim2.csv")
b  <- r[r$sim == "sim2" & r$z1 == "unif_m3_3", ]
shift <- -(b$tau_mix - b$tau_super)   # sign-aligned with the paper (see decisions D5)
difb  <- -b$bias_super
out <- do.call(rbind, lapply(c("cm_bias", "gps_bias", "dif_bias"), function(col) {
  y <- t3[[col]]
  data.frame(paper_column = col,
             cor_with_mix_shift = cor(y, shift), sign_agree_mix_shift = sum(sign(y) == sign(shift)),
             cor_with_pop_DIF_bias = cor(y, difb), sign_agree_pop_DIF_bias = sum(sign(y) == sign(difb)))
}))
write.csv(out, "results/stage1/design_evidence_sim2.csv", row.names = FALSE)
print(out, digits = 3, row.names = FALSE)
