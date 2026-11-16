## Stage 1: how was the paper's mis-normalised probability formula used?
## The paper writes p_t = exp(eta_t) / [1 + sum_{t=1}^A exp(eta_t)] with eta_1 = 0,
## which sums to less than 1.  Variants: (a) normalised softmax (what sample()/rmultinom
## would do); (b) the missing mass 1/D goes to t = 1; (c) it goes to t = A.
## Compared on the paper's DIF bias (Table 1 |bias|; Table 3 signed bias, sign-aligned).
## Run from the repository root:  Rscript scripts/stage1_dgp_variants.R
src <- readLines("scripts/stage1_dgp_checks.R")
eval(parse(text = src[1:(grep("^pop_summary", src) - 1)]))
t1 <- read.csv("data/paper/table1_sim1.csv"); t3 <- read.csv("data/paper/table3_sim2.csv")
probs_v <- function(W, sim, v) {
  Q <- exp(eta_mat(W, sim)); D <- 1 + rowSums(Q); P <- Q / D
  if (v == "softmax")   P <- Q / rowSums(Q)
  if (v == "left_to_1") P[, 1] <- P[, 1] + 1 / D
  if (v == "left_to_A") P[, ncol(P)] <- P[, ncol(P)] + 1 / D
  P
}
run <- function(sim, v, N = 3e6) {
  set.seed(2); W <- draw_W(N, "unif_m3_3"); P <- probs_v(W, sim, v); M <- mu(W, sim)
  A <- ncol(P); pr <- t(combn(A, 2))
  cm <- sapply(1:A, function(t) sum(M[, t] * P[, t]) / sum(P[, t])); sm <- colMeans(M)
  b <- (cm[pr[, 1]] - cm[pr[, 2]]) - (sm[pr[, 1]] - sm[pr[, 2]])
  if (sim == "sim1") data.frame(sim, variant = v, compare = "abs bias vs Table 1",
                                rmsd = sqrt(mean((abs(b) - t1$dif_absbias)^2)), cor = cor(abs(b), t1$dif_absbias))
  else data.frame(sim, variant = v, compare = "-bias vs Table 3 signed bias",
                  rmsd = sqrt(mean((-b - t3$dif_bias)^2)), cor = cor(-b, t3$dif_bias))
}
out <- do.call(rbind, lapply(c("sim1", "sim2"), function(s)
  do.call(rbind, lapply(c("softmax", "left_to_1", "left_to_A"), function(v) run(s, v)))))
write.csv(out, "results/stage1/dgp_assignment_variants.csv", row.names = FALSE)
print(out, digits = 3, row.names = FALSE)
