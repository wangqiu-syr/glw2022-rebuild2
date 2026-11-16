## Summary comparison figures, built from tracked result tables only (no run files).
## From the repository root:  Rscript scripts/comparison_figures.R
## Writes figures/summary_vs_paper.png and figures/sim2_bias_by_treatment_vs_paper.png.
suppressPackageStartupMessages(library(ggplot2))
COL_OURS <- "#2a78d6"; COL_PAPER <- "#eb6834"; SURF <- "#fcfcfb"; INK2 <- "#52514e"; LINK <- "#c9c8c2"
theme_viz <- theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3), panel.grid.minor = element_blank(),
        panel.spacing = unit(1.2, "lines"),
        strip.text = element_text(face = "bold", hjust = 0), legend.position = "top", legend.justification = "left",
        plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = INK2),
        plot.caption = element_text(colour = INK2, hjust = 0))
src_scale <- scale_colour_manual(values = c("This rebuild" = COL_OURS, "Paper" = COL_PAPER), name = NULL)
EST <- c("DIF", "CM", "GPS")

## ---- Figure A: paper vs rebuild at a glance
ag <- read.csv("results/agreement_with_paper.csv")
cs <- read.csv("results/coverage_summary.csv")
SIM <- c(sim1 = "Simulation 1", sim2 = "Simulation 2")
M <- c(b = "Mean |bias| over pairs", m = "Mean MSE over pairs (CM, GPS)", c = "Mean coverage of 95% intervals (%)")
rows <- list()
for (e in EST) {
  a <- ag[ag$estimator == e, ]
  rows[[length(rows) + 1]] <- data.frame(sim = unname(SIM[c("sim1", "sim1", "sim2", "sim2")]), measure = unname(M["b"]), estimator = e,
    source = c("This rebuild", "Paper"), value = c(a$sim1_mean_absbias_ours, a$sim1_mean_absbias_paper,
                                                   a$sim2_mean_absbias_ours, a$sim2_mean_absbias_paper))
  if (e != "DIF") {
    rows[[length(rows) + 1]] <- data.frame(sim = unname(SIM["sim2"]), measure = unname(M["m"]), estimator = e,
      source = c("This rebuild", "Paper"), value = c(a$sim2_mean_mse_ours, a$sim2_mean_mse_paper))
  }
  for (s in names(SIM)) {
    x <- cs[cs$sim == s & cs$estimator == e, ]
    rows[[length(rows) + 1]] <- data.frame(sim = unname(SIM[s]), measure = unname(M["c"]), estimator = e, source = "This rebuild",
                                           value = 100 * x$cov_super_mean)
    if (!is.na(x$paper_mean)) rows[[length(rows) + 1]] <- data.frame(sim = unname(SIM[s]), measure = unname(M["c"]), estimator = e,
                                                                     source = "Paper", value = 100 * x$paper_mean)
  }
}
## Sim 1 MSE: the paper reports none; show ours alone
ov <- read.csv("results/estimator_overview.csv"); ov <- ov[ov$truth == "tau_super" & ov$sim == "sim1", ]
for (e in c("CM", "GPS")) rows[[length(rows) + 1]] <- data.frame(sim = unname(SIM["sim1"]), measure = unname(M["m"]), estimator = e,
                                                                 source = "This rebuild", value = ov$mse[ov$estimator == e])
d <- do.call(rbind, rows); rownames(d) <- NULL
d$measure <- factor(d$measure, M); d$estimator <- factor(d$estimator, rev(EST))
d$source <- factor(d$source, c("This rebuild", "Paper"))
seg <- reshape(d, idvar = c("sim", "measure", "estimator"), timevar = "source", direction = "wide")
seg <- seg[complete.cases(seg), ]
ref <- data.frame(measure = factor(M["c"], M), x = 95)
gA <- ggplot(d, aes(value, estimator)) +
  geom_vline(data = ref, aes(xintercept = x), colour = INK2, linewidth = 0.4, linetype = "dashed") +
  geom_segment(data = seg, aes(x = `value.This rebuild`, xend = value.Paper, y = estimator, yend = estimator),
               colour = LINK, linewidth = 1.6) +
  geom_point(aes(colour = source), size = 3.4) +
  facet_grid(sim ~ measure, scales = "free_x") + src_scale + theme_viz +
  labs(x = NULL, y = NULL, title = "This rebuild vs. the paper at a glance",
       subtitle = "Each row: one estimator; the grey bar joins this rebuild's value to the paper's. Dashed line: nominal 95%.",
       caption = paste("Bias against the population effect, 500 replications per simulation. The paper reports no MSE for Simulation 1",
                       "(only this rebuild is shown) and no coverage for DIF.\nDIF's MSE is omitted: the paper's DIF 'MSE' is not an MSE (decisions D15)."))
ggsave("figures/summary_vs_paper.png", gA, width = 10, height = 5.2, dpi = 200, bg = SURF)

## ---- Figure B: Sim 2 bias by treatment (additive components)
## bias(t,t') = h_t - h_t', so the 36 pairwise biases of each estimator reduce to
## nine per-treatment components; centred to mean zero (only differences matter).
dec <- read.csv("results/sim2/bias_decomposition_vs_paper.csv")
dB <- do.call(rbind, lapply(EST, function(e) {
  o <- dec[[paste0(e, "_ours")]]; p <- dec[[paste0(e, "_paper_aligned")]]
  rbind(data.frame(estimator = e, t = dec$t, source = "This rebuild", value = o - mean(o)),
        data.frame(estimator = e, t = dec$t, source = "Paper", value = p - mean(p)))
}))
dB$estimator <- factor(dB$estimator, EST); dB$source <- factor(dB$source, c("This rebuild", "Paper"))
gB <- ggplot(dB, aes(t, value, colour = source)) +
  geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
  geom_line(linewidth = 0.7) + geom_point(size = 2.4) +
  facet_wrap(~estimator, nrow = 1, scales = "free_y") + src_scale + theme_viz +
  scale_x_continuous(breaks = 1:9, minor_breaks = NULL) +
  labs(x = "treatment combination t", y = "bias component h_t (centred)",
       title = "Simulation 2: where the bias differences sit, by treatment combination",
       subtitle = "Every pairwise bias is h_t - h_t'. Where the two lines diverge, that treatment drives the differences.",
       caption = "Components fitted to the 36 pairwise biases (exact to rounding for the paper's Table 3). Paper's signs multiplied by -1 (decisions D5).")
ggsave("figures/sim2_bias_by_treatment_vs_paper.png", gB, width = 10, height = 4.2, dpi = 200, bg = SURF)
cat("wrote figures/summary_vs_paper.png and figures/sim2_bias_by_treatment_vs_paper.png\n")
