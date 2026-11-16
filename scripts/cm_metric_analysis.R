## CM distance-metric sensitivity (D9): bias and MSE of CM under three metrics,
## against the paper. From the repository root, after
##   Rscript scripts/run_extra.R cm_metrics sim1,sim2 1 500 results/runs_cm_metrics 4
## run:  Rscript scripts/cm_metric_analysis.R
source("R/load.R")
suppressPackageStartupMessages(library(ggplot2))
out <- "results/cm_metric_sensitivity"; dir.create(out, showWarnings = FALSE)
LAB <- c(CM = "blocks: Mahalanobis(X) + Mahalanobis(Z) [primary]",
         CM_stacked = "stacked: Mahalanobis(X, Z)",
         CM_xmahal_zeuclid = "Mahalanobis(X) + Euclidean(Z)")
paper1 <- read.csv("data/paper/table1_sim1.csv"); paper3 <- read.csv("data/paper/table3_sim2.csv")
f4 <- function(x) formatC(x, format = "f", digits = 4)
md_table <- function(df) c(paste0("| ", paste(names(df), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
  apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
perf <- list(); check <- c()
for (s in names(SIMS)) {
  ex <- collect_extra(file.path("results/runs_cm_metrics", s))
  st3 <- read.csv(file.path("results", s, "estimates.csv"))
  a <- merge(ex$estimates[ex$estimates$estimator == "CM", ], st3[st3$estimator == "CM", ],
             by = c("rep", "t", "tprime"))
  check[s] <- max(abs(a$est.x - a$est.y))
  tr <- read.csv(sprintf("results/truth/%s_truth.csv", s))
  perf[[s]] <- rbind(performance(ex$estimates, tr, "tau_super"),
                     performance(ex$estimates, tr, "tau_mix"))
  write.csv(perf[[s]], file.path(out, sprintf("performance_%s.csv", s)), row.names = FALSE)
}
summ <- do.call(rbind, lapply(names(CM_METRICS), function(e) {
  p1 <- perf$sim1[perf$sim1$estimator == e & perf$sim1$truth == "tau_super", ]
  p2 <- perf$sim2[perf$sim2$estimator == e & perf$sim2$truth == "tau_super", ]
  m2 <- perf$sim2[perf$sim2$estimator == e & perf$sim2$truth == "tau_mix", ]
  q2 <- -paper3$cm_bias
  data.frame(metric = e,
             sim1_mean_absbias = mean(p1$abs_bias), sim1_paper = mean(paper1$cm_absbias),
             sim1_cor_absbias = cor(p1$abs_bias, paper1$cm_absbias), sim1_mean_mse = mean(p1$mse),
             sim2_mean_absbias = mean(p2$abs_bias), sim2_paper = mean(abs(q2)),
             sim2_cor_bias = cor(p2$bias, q2), sim2_slope_paper_on_ours = unname(coef(lm(q2 ~ p2$bias))[2]),
             sim2_rmsd = sqrt(mean((p2$bias - q2)^2)), sim2_mean_mse = mean(p2$mse),
             sim2_paper_mse = mean(paper3$cm_mse), sim2_cor_mse = cor(p2$mse, paper3$cm_mse),
             sim2_mean_absbias_vs_mix = mean(m2$abs_bias))
}))
write.csv(summ, file.path(out, "summary.csv"), row.names = FALSE)

## figure: Sim 2 bias by pair, one panel per metric
d <- do.call(rbind, lapply(names(CM_METRICS), function(e) {
  p <- perf$sim2[perf$sim2$estimator == e & perf$sim2$truth == "tau_super", ]
  rbind(data.frame(metric = unname(LAB[e]), idx = 1:36, source = "This rebuild", value = p$bias, se = p$mcse_bias),
        data.frame(metric = unname(LAB[e]), idx = 1:36, source = "Paper CM (sign-aligned)", value = -paper3$cm_bias, se = NA))
}))
d$metric <- factor(d$metric, LAB)
g <- ggplot(d, aes(idx, value, colour = source)) +
  geom_hline(yintercept = 0, colour = "#52514e", linewidth = 0.3) +
  geom_ribbon(data = d[d$source == "This rebuild", ], aes(ymin = value - 1.96 * se, ymax = value + 1.96 * se),
              fill = "#2a78d6", alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  facet_wrap(~metric, ncol = 1) +
  scale_colour_manual(values = c("This rebuild" = "#2a78d6", "Paper CM (sign-aligned)" = "#eb6834"), name = NULL) +
  scale_x_continuous(breaks = c(1, seq(4, 36, 4)), minor_breaks = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "#fcfcfb", colour = NA),
        panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3), panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0), legend.position = "top", legend.justification = "left",
        plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = "#52514e")) +
  labs(x = "pair index (Table 3 order)", y = "bias against tau_super (= 0)",
       title = "Simulation 2: CM bias under three distance metrics, vs. the paper",
       subtitle = "Shaded band: 95% Monte Carlo interval for this rebuild; 500 replications on the Stage-3 data sets")
ggsave("figures/cm_metric_sensitivity_sim2.png", g, width = 9.5, height = 7.5, dpi = 200, bg = "#fcfcfb")

L <- c("# CM distance-metric sensitivity (D9)", "",
       "Same 500 + 500 data sets as Stage 3 (same random streams). The primary CM estimates reproduce Stage 3 exactly",
       sprintf("(max |difference|: Sim 1 %.1e, Sim 2 %.1e).", check["sim1"], check["sim2"]), "",
       "Bias against tau_super; paper's Sim 2 signs multiplied by -1 (D5). Paper CM: Sim 1 mean |bias| 0.1432; Sim 2 mean |bias| 0.1394, mean MSE 0.3677.", "",
       md_table(data.frame(metric = unname(LAB[summ$metric]),
                           `S1 mean abs bias` = f4(summ$sim1_mean_absbias), `S1 cor with paper` = f4(summ$sim1_cor_absbias),
                           `S1 mean MSE` = f4(summ$sim1_mean_mse),
                           `S2 mean abs bias` = f4(summ$sim2_mean_absbias), `S2 cor with paper` = f4(summ$sim2_cor_bias),
                           `S2 slope paper~ours` = f4(summ$sim2_slope_paper_on_ours),
                           `S2 mean MSE` = f4(summ$sim2_mean_mse), `S2 mean abs bias vs tau_mix` = f4(summ$sim2_mean_absbias_vs_mix),
                           check.names = FALSE)), "",
       "Per-pair values: `performance_sim1.csv`, `performance_sim2.csv`. Figure: `figures/cm_metric_sensitivity_sim2.png`.")
writeLines(L, file.path(out, "report.md")); cat(L, sep = "\n")
