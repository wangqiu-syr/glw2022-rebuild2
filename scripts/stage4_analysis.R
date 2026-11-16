## Stage 4: interval coverage (D26). From the repository root, after
##   Rscript scripts/run_extra.R stage4 sim1,sim2 1 500 results/runs_stage4 4
## run:  Rscript scripts/stage4_analysis.R
source("R/load.R")
suppressPackageStartupMessages(library(ggplot2))
Z <- qnorm(0.975)
paper1 <- read.csv("data/paper/table1_sim1.csv"); paper3 <- read.csv("data/paper/table3_sim2.csv")
f1 <- function(x) formatC(100 * x, format = "f", digits = 1)
f3 <- function(x) formatC(x, format = "f", digits = 3)
md_table <- function(df) c(paste0("| ", paste(names(df), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
  apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))

cov_all <- list(); check <- c(); conv <- list()
for (s in names(SIMS)) {
  ex <- collect_extra(file.path("results/runs_stage4", s))
  stopifnot(length(unique(ex$diag$rep)) == N_REPS)
  st3 <- read.csv(file.path("results", s, "estimates.csv"))
  a <- merge(ex$estimates, st3, by = c("rep", "estimator", "t", "tprime"))
  check[s] <- max(abs(a$est.x - a$est.y))
  tr <- read.csv(sprintf("results/truth/%s_truth.csv", s))
  d <- merge(merge(ex$estimates, ex$var[, c("rep", "estimator", "t", "tprime", "var_noise", "var_sate", "var_full")],
                   by = c("rep", "estimator", "t", "tprime")),
             ex$sate[, c("rep", "t", "tprime", "sate")], by = c("rep", "t", "tprime"))
  d <- merge(d, tr[, c("t", "tprime", "tau_super", "tau_mix")], by = c("t", "tprime"))
  write.csv(d[order(d$estimator, d$t, d$tprime, d$rep), ], file.path("results", s, "intervals.csv"), row.names = FALSE)
  key <- interaction(d$estimator, d$t, d$tprime, drop = TRUE)
  cv <- do.call(rbind, lapply(split(d, key), function(g) {
    se <- sqrt(g$var_full); sd_mc <- sd(g$est)
    data.frame(sim = s, estimator = g$estimator[1], t = g$t[1], tprime = g$tprime[1],
               cov_super = mean(abs(g$est - g$tau_super) <= Z * se),
               cov_mix = mean(abs(g$est - g$tau_mix) <= Z * se),
               cov_sate_cond = if (all(is.na(g$var_noise))) NA else mean(abs(g$est - g$sate) <= Z * sqrt(g$var_noise)),
               cov_oracle_super = mean(abs(g$est - g$tau_super) <= Z * sd_mc),
               mean_se = mean(se), sd_mc = sd_mc, se_ratio = mean(se) / sd_mc,
               share_var_units = if (all(is.na(g$var_sate))) NA else mean(g$var_sate / g$var_full))
  }))
  cv$estimator <- factor(cv$estimator, ESTIMATORS); cv <- cv[order(cv$estimator, cv$t, cv$tprime), ]
  cv$estimator <- as.character(cv$estimator); rownames(cv) <- NULL
  if (s == "sim1") {
    cv$paper_coverage <- NA
    for (e in c("CM", "GPS")) cv$paper_coverage[cv$estimator == e] <- paper1[[paste0(tolower(e), "_coverage")]] / 100
  } else {
    cv$paper_coverage <- NA
    for (e in c("CM", "GPS")) cv$paper_coverage[cv$estimator == e] <- paper3[[paste0(tolower(e), "_coverage")]] / 100
  }
  write.csv(cv, file.path("results", s, "coverage.csv"), row.names = FALSE)
  cov_all[[s]] <- cv
}
cv <- do.call(rbind, cov_all)

## summary by estimator
sm <- do.call(rbind, lapply(split(cv, interaction(cv$sim, cv$estimator, drop = TRUE)), function(g) data.frame(
  sim = g$sim[1], estimator = g$estimator[1], pairs = nrow(g),
  cov_super_mean = mean(g$cov_super), cov_super_min = min(g$cov_super), cov_super_max = max(g$cov_super),
  cov_mix_mean = mean(g$cov_mix), cov_sate_mean = mean(g$cov_sate_cond), cov_oracle_mean = mean(g$cov_oracle_super),
  se_ratio_mean = mean(g$se_ratio), se_ratio_min = min(g$se_ratio), se_ratio_max = max(g$se_ratio),
  paper_mean = mean(g$paper_coverage), paper_min = min(g$paper_coverage), paper_max = max(g$paper_coverage))))
sm$estimator <- factor(sm$estimator, ESTIMATORS); sm <- sm[order(sm$sim, sm$estimator), ]
sm$estimator <- as.character(sm$estimator); rownames(sm) <- NULL
write.csv(sm, "results/coverage_summary.csv", row.names = FALSE)

## ---------------------------------------------------------------- figures
COL_OURS <- "#2a78d6"; COL_PAPER <- "#eb6834"; SURF <- "#fcfcfb"; INK2 <- "#52514e"
theme_viz <- theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3), panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0), legend.position = "top", legend.justification = "left",
        plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = INK2),
        plot.caption = element_text(colour = INK2, hjust = 0))
band <- Z * sqrt(0.95 * 0.05 / N_REPS)
dd <- do.call(rbind, lapply(names(SIMS), function(s) {
  g <- cov_all[[s]][cov_all[[s]]$estimator %in% c("CM", "GPS"), ]
  rbind(data.frame(sim = s, estimator = g$estimator, idx = ave(g$t, g$estimator, FUN = seq_along),
                   source = "This rebuild", value = 100 * g$cov_super),
        data.frame(sim = s, estimator = g$estimator, idx = ave(g$t, g$estimator, FUN = seq_along),
                   source = "Paper", value = 100 * g$paper_coverage))
}))
dd$panel <- factor(paste0(ifelse(dd$sim == "sim1", "Simulation 1", "Simulation 2"), ": ", dd$estimator),
                   c("Simulation 1: CM", "Simulation 1: GPS", "Simulation 2: CM", "Simulation 2: GPS"))
gc <- ggplot(dd, aes(idx, value, colour = source)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 100 * (0.95 - band), ymax = 100 * (0.95 + band),
           fill = "#e6e5e0", alpha = 0.6) +
  geom_hline(yintercept = 95, colour = INK2, linewidth = 0.3) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  facet_wrap(~panel, ncol = 2, scales = "free_x") +
  scale_colour_manual(values = c("This rebuild" = COL_OURS, "Paper" = COL_PAPER), name = NULL) +
  scale_x_continuous(breaks = function(l) if (l[2] < 10) 1:6 else c(1, seq(4, 36, 4)), minor_breaks = NULL) +
  theme_viz +
  labs(x = "pair index (table order)", y = "coverage of nominal 95% interval (%)",
       title = "Coverage of 95% intervals for tau_super, this rebuild vs. the paper",
       subtitle = "This rebuild: analytic variance (D26). Grey band: 95% +/- Monte Carlo error for 500 replications.",
       caption = "The paper does not say how its intervals were built.")
ggsave("figures/coverage_vs_paper.png", gc, width = 9.5, height = 6.6, dpi = 200, bg = SURF)

ds <- cv; ds$estimator <- factor(ds$estimator, rev(ESTIMATORS))
ds$sim <- factor(ds$sim, c("sim1", "sim2"), c("Simulation 1", "Simulation 2"))
ms <- aggregate(se_ratio ~ sim + estimator, ds, mean)
gs <- ggplot(ds, aes(se_ratio, estimator)) +
  geom_vline(xintercept = 1, colour = INK2, linewidth = 0.4, linetype = "dashed") +
  geom_point(colour = "#b9b8b2", size = 1.4, position = position_jitter(height = 0.12, seed = 1)) +
  geom_point(data = ms, colour = COL_OURS, size = 3.2) +
  facet_wrap(~sim, nrow = 1) + theme_viz +
  labs(x = "mean estimated SE / Monte Carlo SD of the estimates", y = NULL,
       title = "Calibration of the estimated standard errors",
       subtitle = "Grey: one point per pairwise effect. Blue: mean over pairs. 1 = calibrated.")
ggsave("figures/se_calibration.png", gs, width = 9.5, height = 4.2, dpi = 200, bg = SURF)

## ---------------------------------------------------------------- report
L <- c("# Stage 4 results: intervals and coverage", "",
       sprintf("Generated %s with R %s. Same 500 + 500 data sets as Stage 3; the estimates reproduce Stage 3 (max |difference| Sim 1 %.1e, Sim 2 %.1e).",
               format(Sys.Date()), getRversion(), check["sim1"], check["sim2"]),
       "Interval: estimate +/- 1.96 x estimated SE, with the analytic variance of decisions D26 (DIF: s_t^2/n_t + s_t'^2/n_t'). Monte Carlo SE of a coverage near 95% is about 1 percentage point.", "",
       "Columns: coverage (%) against tau_super (the paper's truth) as mean [min, max] over pairs; against tau_mix (the estimand under the equal-group design); against the sample-average effect with the conditional variance; the oracle interval (true Monte Carlo SD) against tau_super; mean estimated SE / Monte Carlo SD; the paper's coverage as mean [min, max].", "",
       md_table(data.frame(sim = sm$sim, estimator = sm$estimator,
                           `vs tau_super` = paste0(f1(sm$cov_super_mean), " [", f1(sm$cov_super_min), ", ", f1(sm$cov_super_max), "]"),
                           `vs tau_mix` = f1(sm$cov_mix_mean),
                           `vs sample effect (cond.)` = ifelse(is.na(sm$cov_sate_mean), "-", f1(sm$cov_sate_mean)),
                           `oracle vs tau_super` = f1(sm$cov_oracle_mean),
                           `SE / MC SD` = paste0(f3(sm$se_ratio_mean), " [", f3(sm$se_ratio_min), ", ", f3(sm$se_ratio_max), "]"),
                           paper = ifelse(is.na(sm$paper_mean), "-",
                                          paste0(f1(sm$paper_mean), " [", f1(sm$paper_min), ", ", f1(sm$paper_max), "]")),
                           check.names = FALSE)), "")
for (s in names(SIMS)) {
  g <- cov_all[[s]]; pr <- unique(g[, c("t", "tprime")])
  tab <- data.frame(pair = sprintf("tau(%d,%d)", pr$t, pr$tprime))
  for (e in ESTIMATORS) tab[[e]] <- f1(g$cov_super[g$estimator == e])
  tab[["CM paper"]] <- f1(g$paper_coverage[g$estimator == "CM"])
  tab[["GPS paper"]] <- f1(g$paper_coverage[g$estimator == "GPS"])
  L <- c(L, sprintf("## %s: coverage (%%) against tau_super, by pair", ifelse(s == "sim1", "Simulation 1", "Simulation 2")), "",
         md_table(tab), "")
}
L <- c(L, "Per-pair details (all coverage types, SE ratio, share of variance from the units term): `results/<sim>/coverage.csv`; per-replication intervals: `results/<sim>/intervals.csv`. Figures: `figures/coverage_vs_paper.png`, `figures/se_calibration.png`.")
writeLines(L, "results/stage4_report.md")
cat(L[1:20], sep = "\n")
