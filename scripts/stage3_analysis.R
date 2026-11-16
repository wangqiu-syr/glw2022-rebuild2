## Stage 3: performance (bias, |bias|, MSE) of all estimators at R = 500,
## side by side with the paper's Tables 1 and 3, plus figures.
## From the repository root, after
##   Rscript scripts/run_simulations.R sim1,sim2 1 500 results/runs 4
## run:  Rscript scripts/stage3_analysis.R
source("R/load.R")
suppressPackageStartupMessages(library(ggplot2))

paper1 <- read.csv("data/paper/table1_sim1.csv")
paper3 <- read.csv("data/paper/table3_sim2.csv")
PAPER_EST <- c("DIF", "CM", "GPS")
f4 <- function(x) formatC(x, format = "f", digits = 4)
f3 <- function(x) formatC(x, format = "f", digits = 3)
md_table <- function(df) {
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
## additive decomposition bias(t,t') = h_t - h_t' (h_1 = 0), see spec 3
additive_h <- function(t, tp, y, A) {
  X <- matrix(0, length(y), A); X[cbind(seq_along(y), t)] <- 1; X[cbind(seq_along(y), tp)] <- -1
  c(0, lm.fit(X[, -1, drop = FALSE], y)$coefficients)
}

res <- list()
for (s in names(SIMS)) {
  dir.create(file.path("results", s), showWarnings = FALSE)
  col <- collect_simulation(file.path("results/runs", s))
  stopifnot(length(unique(col$diag$rep)) == N_REPS)
  tr <- read.csv(sprintf("results/truth/%s_truth.csv", s))
  write.csv(col$estimates, file.path("results", s, "estimates.csv"), row.names = FALSE)
  write.csv(col$diag, file.path("results", s, "replication_diagnostics.csv"), row.names = FALSE)
  write.csv(col$sate, file.path("results", s, "sample_average_effects.csv"), row.names = FALSE)
  perf <- rbind(performance(col$estimates, tr, "tau_super"), performance(col$estimates, tr, "tau_mix"))
  write.csv(perf, file.path("results", s, "performance.csv"), row.names = FALSE)
  vs <- aggregate(varsum_sample ~ t + tprime, col$varsum, mean)
  res[[s]] <- list(perf = perf, truth = tr, diag = col$diag, varsum = vs)
}

## ---------------------------------------------------------------- Sim 1
p1 <- res$sim1$perf[res$sim1$perf$truth == "tau_super", ]
get <- function(p, e, what) p[p$estimator == e, what]
cmp1 <- data.frame(effect = paper1$effect, t = paper1$t, tprime = paper1$tprime,
                   tau = get(p1, "DIF", "tau"))
for (e in ESTIMATORS) {
  cmp1[[paste0(e, "_bias")]] <- get(p1, e, "bias")
  cmp1[[paste0(e, "_absbias")]] <- get(p1, e, "abs_bias")
  cmp1[[paste0(e, "_mse")]] <- get(p1, e, "mse")
  cmp1[[paste0(e, "_mcse_bias")]] <- get(p1, e, "mcse_bias")
}
for (e in PAPER_EST) cmp1[[paste0("paper_", e, "_absbias")]] <- paper1[[paste0(tolower(e), "_absbias")]]
cmp1$paper_CM_coverage <- paper1$cm_coverage; cmp1$paper_GPS_coverage <- paper1$gps_coverage
write.csv(cmp1, "results/sim1/comparison_with_paper.csv", row.names = FALSE)

## ---------------------------------------------------------------- Sim 2
p2 <- res$sim2$perf[res$sim2$perf$truth == "tau_super", ]
cmp2 <- data.frame(effect = paper3$effect, t = paper3$t, tprime = paper3$tprime, idx = seq_len(36))
for (e in ESTIMATORS) {
  cmp2[[paste0(e, "_bias")]] <- get(p2, e, "bias")
  cmp2[[paste0(e, "_mse")]] <- get(p2, e, "mse")
  cmp2[[paste0(e, "_mcse_bias")]] <- get(p2, e, "mcse_bias")
}
for (e in PAPER_EST) {
  cmp2[[paste0("paper_", e, "_bias")]] <- paper3[[paste0(tolower(e), "_bias")]]
  cmp2[[paste0("paper_", e, "_bias_aligned")]] <- -paper3[[paste0(tolower(e), "_bias")]]   # D5
  cmp2[[paste0("paper_", e, "_mse")]] <- paper3[[paste0(tolower(e), "_mse")]]
}
cmp2$DIF_varsum_sample <- res$sim2$varsum$varsum_sample[match(paste(cmp2$t, cmp2$tprime),
                                                         paste(res$sim2$varsum$t, res$sim2$varsum$tprime))]
cmp2$DIF_varsum_pop <- res$sim2$truth$varsum
cmp2$DIF_bias_pop <- res$sim2$truth$dif_bias_pop
cmp2$paper_CM_coverage <- paper3$cm_coverage; cmp2$paper_GPS_coverage <- paper3$gps_coverage
write.csv(cmp2, "results/sim2/comparison_with_paper.csv", row.names = FALSE)

## per-treatment decomposition of bias (ours vs paper, sign-aligned)
dec <- data.frame(t = 1:9)
for (e in PAPER_EST) {
  dec[[paste0(e, "_ours")]] <- additive_h(cmp2$t, cmp2$tprime, cmp2[[paste0(e, "_bias")]], 9)
  dec[[paste0(e, "_paper_aligned")]] <- additive_h(cmp2$t, cmp2$tprime, cmp2[[paste0("paper_", e, "_bias_aligned")]], 9)
}
write.csv(dec, "results/sim2/bias_decomposition_vs_paper.csv", row.names = FALSE)

## ---------------------------------------------------------------- agreement summary
agree <- do.call(rbind, lapply(PAPER_EST, function(e) {
  o1 <- cmp1[[paste0(e, "_absbias")]]; q1 <- cmp1[[paste0("paper_", e, "_absbias")]]
  se1 <- cmp1[[paste0(e, "_mcse_bias")]]
  o2 <- cmp2[[paste0(e, "_bias")]]; q2 <- cmp2[[paste0("paper_", e, "_bias_aligned")]]
  se2 <- cmp2[[paste0(e, "_mcse_bias")]]
  m2 <- cmp2[[paste0(e, "_mse")]]; qm2 <- cmp2[[paste0("paper_", e, "_mse")]]
  ## |ours - paper| within 2 SE of a difference of two MC estimates (paper's MC SE ~ ours)
  data.frame(estimator = e,
             sim1_mean_absbias_ours = mean(o1), sim1_mean_absbias_paper = mean(q1),
             sim1_cor_absbias = cor(o1, q1), sim1_within_2se = sum(abs(o1 - q1) < 2 * sqrt(2) * se1),
             sim2_mean_absbias_ours = mean(abs(o2)), sim2_mean_absbias_paper = mean(abs(q2)),
             sim2_cor_bias_aligned = cor(o2, q2), sim2_rmsd_bias = sqrt(mean((o2 - q2)^2)),
             sim2_sign_agree = sum(sign(o2) == sign(q2)),
             sim2_within_2se = sum(abs(o2 - q2) < 2 * sqrt(2) * se2),
             sim2_mean_mse_ours = mean(m2), sim2_mean_mse_paper = mean(qm2),
             sim2_median_mse_ratio_paper_over_ours = median(qm2 / m2))
}))
write.csv(agree, "results/agreement_with_paper.csv", row.names = FALSE)

## estimator overview incl. sensitivity variants, both truths
ov <- do.call(rbind, lapply(names(SIMS), function(s) {
  p <- res[[s]]$perf
  a <- aggregate(cbind(abs_bias, mse, sd, mcse_bias) ~ estimator + truth, p, mean)
  mx <- aggregate(abs_bias ~ estimator + truth, p, max); names(mx)[3] <- "max_abs_bias"
  cbind(sim = s, merge(a, mx))
}))
ov$estimator <- factor(ov$estimator, ESTIMATORS)
ov <- ov[order(ov$sim, ov$truth != "tau_super", ov$estimator), ]
write.csv(ov, "results/estimator_overview.csv", row.names = FALSE)

## ---------------------------------------------------------------- figures
COL_OURS <- "#2a78d6"; COL_PAPER <- "#eb6834"; COL_THIRD <- "#1baf7a"
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e0"
theme_viz <- function() theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.background = element_rect(fill = SURF, colour = NA),
        panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
        panel.grid.minor = element_blank(),
        text = element_text(colour = INK), axis.text = element_text(colour = INK2),
        strip.text = element_text(face = "bold", hjust = 0, colour = INK),
        legend.position = "top", legend.justification = "left",
        plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = INK2),
        plot.caption = element_text(colour = INK2, hjust = 0))
src_scale <- scale_colour_manual(values = c("This rebuild" = COL_OURS, "Paper" = COL_PAPER,
                                            "Paper (sign-aligned)" = COL_PAPER,
                                            "This rebuild: Var(Y|t) + Var(Y|t')" = COL_THIRD), name = NULL)
dir.create("figures", showWarnings = FALSE)

## Fig 1: Sim 1 absolute bias
d1 <- do.call(rbind, lapply(PAPER_EST, function(e) rbind(
  data.frame(estimator = e, pair = cmp1$effect, source = "This rebuild",
             value = cmp1[[paste0(e, "_absbias")]], se = cmp1[[paste0(e, "_mcse_bias")]]),
  data.frame(estimator = e, pair = cmp1$effect, source = "Paper",
             value = cmp1[[paste0("paper_", e, "_absbias")]], se = NA))))
d1$estimator <- factor(d1$estimator, PAPER_EST); d1$pair <- factor(d1$pair, cmp1$effect)
g1 <- ggplot(d1, aes(pair, value, colour = source)) +
  geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(value - 1.96 * se, 0), ymax = value + 1.96 * se), width = 0.15,
                linewidth = 0.5, position = position_dodge(0.45), na.rm = TRUE) +
  geom_point(size = 2.6, position = position_dodge(0.45)) +
  facet_wrap(~estimator, nrow = 1) + src_scale + theme_viz() +
  labs(x = NULL, y = "|bias| against tau_super",
       title = "Simulation 1: absolute bias, this rebuild vs. paper's Table 1",
       subtitle = "500 replications each; bars are 95% Monte Carlo intervals for this rebuild",
       caption = "Paper values from Table 1. tau_super is the population effect (docs/spec.md section 4).") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/sim1_absbias_vs_paper.png", g1, width = 9.5, height = 4.2, dpi = 200, bg = SURF)

## Fig 2: Sim 2 signed bias by pair (cf. the paper's Fig. 1)
d2 <- do.call(rbind, lapply(PAPER_EST, function(e) rbind(
  data.frame(estimator = e, idx = cmp2$idx, source = "This rebuild", value = cmp2[[paste0(e, "_bias")]],
             se = cmp2[[paste0(e, "_mcse_bias")]]),
  data.frame(estimator = e, idx = cmp2$idx, source = "Paper (sign-aligned)",
             value = cmp2[[paste0("paper_", e, "_bias_aligned")]], se = NA))))
d2$estimator <- factor(d2$estimator, PAPER_EST)
g2 <- ggplot(d2, aes(idx, value, colour = source)) +
  geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
  geom_ribbon(data = d2[d2$source == "This rebuild", ],
              aes(ymin = value - 1.96 * se, ymax = value + 1.96 * se), fill = COL_OURS, alpha = 0.15,
              colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  facet_wrap(~estimator, ncol = 1, scales = "free_y") + src_scale + theme_viz() +
  scale_x_continuous(breaks = c(1, seq(4, 36, 4)), minor_breaks = NULL) +
  labs(x = "pair index (1 = tau(1,2), ..., 36 = tau(8,9), in Table 3 order)", y = "bias against tau_super (= 0)",
       title = "Simulation 2: bias for the 36 pairwise effects, this rebuild vs. paper's Table 3",
       subtitle = "Paper's signs multiplied by -1 (decisions D5); shaded band: 95% Monte Carlo interval for this rebuild",
       caption = "500 replications. Panels have separate y-scales.")
ggsave("figures/sim2_bias_vs_paper.png", g2, width = 9.5, height = 8, dpi = 200, bg = SURF)

## Fig 3: Sim 2 MSE by pair (log scale)
d3 <- do.call(rbind, lapply(PAPER_EST, function(e) rbind(
  data.frame(estimator = e, idx = cmp2$idx, source = "This rebuild", value = cmp2[[paste0(e, "_mse")]]),
  data.frame(estimator = e, idx = cmp2$idx, source = "Paper", value = cmp2[[paste0("paper_", e, "_mse")]]),
  if (e == "DIF") data.frame(estimator = e, idx = cmp2$idx, source = "This rebuild: Var(Y|t) + Var(Y|t')",
                             value = cmp2$DIF_varsum_sample))))
d3$estimator <- factor(d3$estimator, PAPER_EST)
g3 <- ggplot(d3, aes(idx, value, colour = source)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  facet_wrap(~estimator, ncol = 1, scales = "free_y") + src_scale + theme_viz() +
  scale_y_log10() + scale_x_continuous(breaks = c(1, seq(4, 36, 4)), minor_breaks = NULL) +
  labs(x = "pair index (Table 3 order)", y = "MSE (log scale)",
       title = "Simulation 2: MSE, this rebuild vs. paper's Table 3",
       subtitle = "DIF panel also shows the sum of within-group outcome variances, which the paper's DIF 'MSE' tracks (D15)",
       caption = "500 replications. Panels have separate y-scales.")
ggsave("figures/sim2_mse_vs_paper.png", g3, width = 9.5, height = 8, dpi = 200, bg = SURF)

## Fig 4: scatter, ours vs paper
d4 <- rbind(
  do.call(rbind, lapply(PAPER_EST, function(e) data.frame(
    panel = paste("Sim 1 |bias|:", e), ours = cmp1[[paste0(e, "_absbias")]],
    paper = cmp1[[paste0("paper_", e, "_absbias")]]))),
  do.call(rbind, lapply(PAPER_EST, function(e) data.frame(
    panel = paste("Sim 2 bias (sign-aligned):", e), ours = cmp2[[paste0(e, "_bias")]],
    paper = cmp2[[paste0("paper_", e, "_bias_aligned")]]))))
d4$panel <- factor(d4$panel, unique(d4$panel))
lims <- do.call(rbind, lapply(split(d4, d4$panel), function(g) {
  r <- range(c(g$ours, g$paper)); data.frame(panel = g$panel[1], ours = r, paper = r) }))
g4 <- ggplot(d4, aes(paper, ours)) +
  geom_blank(data = lims) +
  geom_abline(slope = 1, intercept = 0, colour = INK2, linewidth = 0.4, linetype = "dashed") +
  geom_point(colour = COL_OURS, size = 2.2) +
  facet_wrap(~panel, nrow = 2, scales = "free") + theme_viz() +
  labs(x = "paper", y = "this rebuild", title = "Bias: this rebuild against the paper, by pair",
       subtitle = "Dashed line: equality. Each point is one pairwise effect.")
ggsave("figures/bias_scatter_vs_paper.png", g4, width = 9.5, height = 6.4, dpi = 200, bg = SURF)

## Fig 5: estimator overview incl. sensitivity variants
d5 <- do.call(rbind, lapply(names(SIMS), function(s) {
  p <- res[[s]]$perf[res[[s]]$perf$truth == "tau_super", ]
  rbind(data.frame(sim = s, estimator = p$estimator, measure = "|bias| against tau_super", value = p$abs_bias),
        data.frame(sim = s, estimator = p$estimator, measure = "MSE", value = p$mse))
}))
d5$estimator <- factor(d5$estimator, rev(ESTIMATORS)); d5$sim <- factor(d5$sim, c("sim1", "sim2"), c("Simulation 1", "Simulation 2"))
m5 <- aggregate(value ~ sim + estimator + measure, d5, mean)
g5 <- ggplot(d5, aes(value, estimator)) +
  geom_point(colour = "#b9b8b2", size = 1.4, position = position_jitter(height = 0.12, seed = 1)) +
  geom_point(data = m5, colour = COL_OURS, size = 3.2) +
  facet_grid(sim ~ measure, scales = "free_x") + theme_viz() +
  labs(x = NULL, y = NULL, title = "All estimators, including the GPS sensitivity variants",
       subtitle = "Grey: one point per pairwise effect. Blue: mean over pairs. 500 replications.",
       caption = "GPS: main-effects logit (primary). GPS_cs: correctly specified logit. GPS_true: true propensity scores.")
ggsave("figures/estimator_overview.png", g5, width = 9.5, height = 5.2, dpi = 200, bg = SURF)

## ---------------------------------------------------------------- report
dg <- do.call(rbind, lapply(names(SIMS), function(s) {
  d <- res[[s]]$diag
  data.frame(sim = s, reps = nrow(d), sec_per_rep = mean(d$seconds),
             main_conv = sum(d$main_converged), cs_conv = sum(d$cs_converged),
             max_iter = paste(max(d$main_iter), "/", max(d$cs_iter)))
}))
L <- c("# Stage 3 results: bias and MSE at R = 500", "",
       sprintf("Generated %s with R %s. Truth for all comparisons: the population effect tau_super (D4); bias against tau_mix is in `results/<sim>/performance.csv`. Bias = E(tau_hat) - tau with tau(t,t') = E[Y(t) - Y(t')]. The paper's signed Table 3 values are shown both as printed and multiplied by -1 (\"aligned\", D5).",
               format(Sys.Date()), getRversion()), "",
       "## Runs", "", md_table(transform(dg, sec_per_rep = f3(sec_per_rep))), "",
       "Columns main_conv and cs_conv count replications in which the main-effects and correctly specified logit fits converged.", "")

## Sim 1 table
t1 <- data.frame(pair = cmp1$effect, tau = cmp1$tau,
                 `DIF ours` = f4(cmp1$DIF_absbias), `DIF paper` = f4(cmp1$paper_DIF_absbias),
                 `CM ours` = f4(cmp1$CM_absbias), `CM paper` = f4(cmp1$paper_CM_absbias),
                 `GPS ours` = f4(cmp1$GPS_absbias), `GPS paper` = f4(cmp1$paper_GPS_absbias), check.names = FALSE)
t1b <- data.frame(pair = cmp1$effect,
                  `DIF bias` = f4(cmp1$DIF_bias), `DIF MSE` = f4(cmp1$DIF_mse),
                  `CM bias` = f4(cmp1$CM_bias), `CM MSE` = f4(cmp1$CM_mse),
                  `GPS bias` = f4(cmp1$GPS_bias), `GPS MSE` = f4(cmp1$GPS_mse),
                  `GPS_cs bias` = f4(cmp1$GPS_cs_bias), `GPS_cs MSE` = f4(cmp1$GPS_cs_mse),
                  `GPS_true bias` = f4(cmp1$GPS_true_bias), `GPS_true MSE` = f4(cmp1$GPS_true_mse), check.names = FALSE)
mc1 <- max(sapply(PAPER_EST, function(e) max(cmp1[[paste0(e, "_mcse_bias")]])))
L <- c(L, "## Simulation 1 (A = 4, n_t = 100)", "",
       "Absolute bias, this rebuild vs. paper Table 1:", "", md_table(t1), "",
       sprintf("Largest Monte Carlo SE of a bias in this table: %.4f.", mc1), "",
       "This rebuild, signed bias and MSE (the paper reports neither for Sim 1):", "", md_table(t1b), "")

## Sim 2 tables
t2 <- data.frame(pair = cmp2$effect,
                 `DIF ours` = f4(cmp2$DIF_bias), `DIF paper (aligned)` = f4(cmp2$paper_DIF_bias_aligned),
                 `CM ours` = f4(cmp2$CM_bias), `CM paper (aligned)` = f4(cmp2$paper_CM_bias_aligned),
                 `GPS ours` = f4(cmp2$GPS_bias), `GPS paper (aligned)` = f4(cmp2$paper_GPS_bias_aligned),
                 check.names = FALSE)
t2m <- data.frame(pair = cmp2$effect,
                  `DIF ours` = f4(cmp2$DIF_mse), `DIF paper` = f4(cmp2$paper_DIF_mse),
                  `DIF Var(Y|t)+Var(Y|t')` = formatC(cmp2$DIF_varsum_sample, format = "f", digits = 1),
                  `CM ours` = f4(cmp2$CM_mse), `CM paper` = f4(cmp2$paper_CM_mse),
                  `GPS ours` = f4(cmp2$GPS_mse), `GPS paper` = f4(cmp2$paper_GPS_mse), check.names = FALSE)
mc2 <- sapply(PAPER_EST, function(e) max(cmp2[[paste0(e, "_mcse_bias")]]))
L <- c(L, "## Simulation 2 (A = 9, n_t = 200; tau_super = 0 for all pairs)", "",
       "Bias, this rebuild vs. paper Table 3 (paper's signs multiplied by -1):", "", md_table(t2), "",
       sprintf("Largest Monte Carlo SE of a bias: DIF %.4f, CM %.4f, GPS %.4f.", mc2[1], mc2[2], mc2[3]), "",
       "MSE, this rebuild vs. paper Table 3. The column Var(Y|t)+Var(Y|t') is the replication mean of the sum of the two within-group sample variances, the quantity the paper's DIF 'MSE' tracks (D15):", "",
       md_table(t2m), "")

## sensitivity / overview
ovt <- ov[ov$truth == "tau_super", ]
L <- c(L, "## All estimators, including the GPS sensitivity variants", "",
       "Means over pairs (against tau_super):", "",
       md_table(data.frame(sim = ovt$sim, estimator = ovt$estimator, `mean |bias|` = f4(ovt$abs_bias),
                           `max |bias|` = f4(ovt$max_abs_bias), `mean MSE` = f4(ovt$mse),
                           `mean SD` = f4(ovt$sd), check.names = FALSE)), "")
ovm <- ov[ov$truth == "tau_mix", ]
L <- c(L, "Against tau_mix, the estimand that CM and GPS target under the equal-group-size design (D4):", "",
       md_table(data.frame(sim = ovm$sim, estimator = ovm$estimator, `mean |bias|` = f4(ovm$abs_bias),
                           `max |bias|` = f4(ovm$max_abs_bias), `mean MSE` = f4(ovm$mse), check.names = FALSE)), "")

## agreement
L <- c(L, "## Agreement with the paper", "",
       md_table(data.frame(estimator = agree$estimator,
                           `S1 mean abs bias ours / paper` = paste(f4(agree$sim1_mean_absbias_ours), "/", f4(agree$sim1_mean_absbias_paper)),
                           `S1 cor` = f3(agree$sim1_cor_absbias), `S1 within 2 SE` = paste0(agree$sim1_within_2se, "/6"),
                           `S2 mean abs bias ours / paper` = paste(f4(agree$sim2_mean_absbias_ours), "/", f4(agree$sim2_mean_absbias_paper)),
                           `S2 cor (aligned)` = f3(agree$sim2_cor_bias_aligned),
                           `S2 signs agree` = paste0(agree$sim2_sign_agree, "/36"),
                           `S2 within 2 SE` = paste0(agree$sim2_within_2se, "/36"),
                           `S2 mean MSE ours / paper` = paste(f4(agree$sim2_mean_mse_ours), "/", f4(agree$sim2_mean_mse_paper)),
                           check.names = FALSE)), "",
       "Columns: S1 = Simulation 1 (absolute biases), S2 = Simulation 2 (signed biases, paper sign-aligned). \"Within 2 SE\": |ours - paper| < 2 x sqrt(2) x our Monte Carlo SE, i.e. roughly two SEs of a difference between two independent 500-replication estimates.", "")
writeLines(L, "results/stage3_report.md")
cat("wrote results/stage3_report.md and figures/\n")
print(agree, digits = 3)
print(ov[ov$truth == "tau_super", ], digits = 3)
