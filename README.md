# glw2022-rebuild2

An independent rebuild, in R, of the two simulation studies in

> Guo, S., Liu, J., & Wang, Q. (2022). Effective learning during COVID-19:
> Multilevel covariates matching and propensity score matching.
> *Annals of Data Science, 9*(5), 967–982.

It was built from the paper alone. No code, replications or discussion of the
paper were used, nor the cited works' implementations. The matching estimators
and the multinomial logit are written from scratch in base R. The paper's PDF is
not in this repository (copyright).

## Status

| Stage | Content | Status |
|---|---|---|
| 1 | Specification, decisions, paper's tables as CSV | done (reviewed) |
| 2 | Code in `R/`, tests, population truths, pilot run | done (reviewed) |
| 3 | Full runs (500 replications each): bias, \|bias\|, MSE against the paper | done (reviewed); CM metric sensitivity added |
| 4 | Intervals and coverage | **done, awaiting review** |

## Layout

| Path | Content |
|---|---|
| [`docs/spec.md`](docs/spec.md) | DGP, estimators, true effects, performance measures, replications |
| [`docs/decisions.md`](docs/decisions.md) | 26 recorded choices where the paper is unclear or inconsistent, or where implementation needed a choice |
| [`data/paper/`](data/paper) | the paper's Tables 1 and 3 as CSV (checked) |
| `R/config.R` | simulation settings (A, n per group, master seeds, estimators) |
| `R/dgp.R` | covariates, treatment model, outcomes, group-filling sampler |
| `R/mlogit.R` | multinomial logit by Newton–Raphson; GPS design matrices |
| `R/estimators.R` | DIF, CM (Mahalanobis 1-NN), GPS matching |
| `R/truth.R` | τ_super (analytic), τ_mix and other population quantities (Monte Carlo integration) |
| `R/rng.R`, `R/run.R` | per-replication random streams; resumable parallel runner; result collection |
| `R/performance.R` | bias, \|bias\|, MSE, SD, Monte Carlo SEs |
| `tests/testthat/` | 111 tests, including hand-worked examples |
| `scripts/` | drivers: tests, truths, simulation runs, pilot summary, stage-1 checks |
| `results/truth/` | population truths per simulation |
| `results/pilot/` | 20-replication pilot: per-replication files, performance preview, timing |
| `results/runs/` | full runs: one file per replication (500 per simulation) |
| `results/sim1/`, `results/sim2/` | per-replication estimates, performance (against both truths), comparison with the paper, diagnostics |
| [`results/stage3_report.md`](results/stage3_report.md) | all Stage 3 tables, side by side with the paper |
| `results/agreement_with_paper.csv`, `results/estimator_overview.csv` | summary tables |
| `figures/` | comparisons with the paper's tables |
| `results/runs_cm_metrics/`, `results/cm_metric_sensitivity/` | CM distance-metric sensitivity (D9) |
| `R/variance.R`, `results/runs_stage4/` | Stage 4 variance estimators (D26) and per-replication runs |
| [`results/stage4_report.md`](results/stage4_report.md), `results/<sim>/coverage.csv`, `results/coverage_summary.csv` | coverage tables; per-replication intervals in `results/<sim>/intervals.csv.gz` |
| `results/stage1/` | stage-1 diagnostic outputs |

Estimators computed in every replication: `DIF`, `CM`, `GPS` (main-effects
logit, primary), plus the sensitivity variants `GPS_cs` (correctly specified
logit) and `GPS_true` (true propensity scores).

## Reproducing

R 4.3.3 from the Ubuntu 24.04 archive (`r-base-core`, `r-cran-testthat`; `r-cran-nnet`
and `r-cran-mass` for one test and the stage-1 scripts). From the repository root:

```sh
Rscript scripts/run_tests.R                                   # 111 tests (~10 s)
Rscript scripts/compute_truth.R                               # population truths (~2.5 min, 4 cores)
Rscript scripts/run_simulations.R sim1,sim2 1 20 results/pilot/runs 4   # pilot (~10 s)
Rscript scripts/pilot_summary.R                               # results/pilot/pilot_summary.md
Rscript scripts/run_simulations.R sim1,sim2 1 500 results/runs 4  # full runs (~2.5 min)
Rscript scripts/stage3_analysis.R                             # tables, report and figures (~1 min)
Rscript scripts/run_extra.R cm_metrics sim1,sim2 1 500 results/runs_cm_metrics 4  # CM metrics (~2.5 min)
Rscript scripts/cm_metric_analysis.R
Rscript scripts/run_extra.R stage4 sim1,sim2 1 500 results/runs_stage4 4      # variances (~3.5 min)
Rscript scripts/stage4_analysis.R                             # coverage tables and figures
```

`run_simulations.R` saves each replication as it finishes and skips finished ones
when rerun, so an interrupted run can be resumed. Replication *r* uses its own
L'Ecuyer-CMRG stream, so it reproduces whether run alone or in a batch.

Stage-1 checks (population-level comparison with the paper's DIF column):

```sh
Rscript scripts/stage1_check_tables.R; Rscript scripts/stage1_dgp_checks.R
Rscript scripts/stage1_dgp_variants.R; Rscript scripts/stage1_design_evidence.R
```

## Main results (Stage 3, 500 replications per simulation)

Bias is $E(\hat\tau)-\tau$ against the population effect $\tau_{super}$ (D4). The
paper's signed biases are compared after multiplying them by −1 (D5). Full
tables are in [`results/stage3_report.md`](results/stage3_report.md).

Means over pairs:

| | Sim 1: mean \|bias\| ours / paper | Sim 1: mean MSE ours | Sim 2: mean \|bias\| ours / paper | Sim 2: mean MSE ours / paper |
|---|---|---|---|---|
| DIF | 0.58 / 0.44 | 1.78 | 0.96 / 1.08 | 2.30 / 182.9 (not an MSE, D15) |
| CM | 0.17 / 0.14 | 0.36 | 0.30 / 0.14 | 0.29 / 0.37 |
| GPS | 0.04 / 0.29 | 0.66 | 0.11 / 0.31 | 0.61 / 1.92 |
| GPS_cs (sensitivity) | 0.04 | 0.71 | 0.10 | 0.84 |
| GPS_true (sensitivity) | 0.11 | 1.58 | 0.18 | 1.66 |

The largest Monte Carlo SE of any bias here is 0.06.

**Reproduced.**
* Both matching estimators remove most of the bias of the sample difference, in
  both simulations.
* CM has a smaller MSE than GPS in both simulations (6 of 6 pairs and 33 of 36
  pairs). This is the paper's "covariates matching is more efficient".
* Sim 2's DIF bias pattern matches the paper closely: correlation 0.96, 30 of 36
  signs agree (`figures/sim2_bias_vs_paper.png`). Our simulated DIF also agrees
  with its population value, so the DGP is implemented as specified.
* CM's MSE is of the same size as the paper's (mean 0.29 vs 0.37).

**Not reproduced.**
* **Which matching estimator is less biased.** In the paper, CM is less biased
  than GPS for 6 of 6 pairs (Sim 1) and 32 of 36 (Sim 2). Here it is the other
  way round: CM is less biased for 0 of 6 and 6 of 36. The paper's GPS bias is about 2.3× ours, and its GPS MSE about 3×
  ours, although the patterns correlate (0.83 in Sim 2). Its CM bias is about a
  third of ours (correlation 0.73).
* **DIF "MSE".** The paper's Table 3 values are 0.81–1.00 × the sum of the two
  within-group outcome variances observed in our runs, i.e. not an MSE
  (`figures/sim2_mse_vs_paper.png`, D15).
* **Some DIF biases.** The pairs flagged at Stage 1 still differ (D20), e.g. Sim 1
  τ(1,3): 1.14 here vs 0.47 in the paper. Because our DIF matches its population
  value, the gap lies in the paper's numbers or in an unstated detail of its DGP.

**CM distance-metric sensitivity** (D9, `results/cm_metric_sensitivity/report.md`).
Rerunning CM with one Mahalanobis distance on all seven covariates, or with
Mahalanobis on X plus Euclidean on Z, changes little. Sim 2 mean |bias| is
0.27–0.30 under all three metrics (paper 0.14), so the metric does not explain
the gap.

**GPS sensitivity variants.**
* The correctly specified logit (`GPS_cs`) gives about the same bias as the
  main-effects logit and a larger MSE (×1.08 in Sim 1, ×1.37 in Sim 2). So in this rebuild,
  misspecification of the propensity model is not what makes GPS less efficient
  than CM, contrary to the paper's explanation.
* Matching on the **true** scores (`GPS_true`) is worse, with larger bias and
  2.4× (Sim 1) to 2.7× (Sim 2) the MSE. This follows from the equal-group-size design (D3). Within
  the sample, the probability of treatment $t$ is
  $q_t(x) \propto p_t(x)/\pi_t$, which is not a function of the population score
  $p_t(x)$ alone, so matching on $p_t$ does not balance the sample. A logit fitted
  to the sample estimates $q_t$, and the correctly specified one is exactly right
  for $q_t$, because only the intercepts shift. Against $\tau_{mix}$, the estimand
  the matching estimators target in this design, `GPS_cs` has the smallest bias
  of all estimators (0.04 and 0.07).
* CM keeps a clear bias even against $\tau_{mix}$ (33 of 36 pairs significant in
  Sim 2). This is consistent with residual covariate differences in one-to-one
  matching on seven covariates with 100–200 donors per group.

## Intervals and coverage (Stage 4)

The paper does not say how its intervals were built. Here the interval is the
estimate ± 1.96 × an analytic SE derived from the estimator's linear-weight form
(D26). That SE has two parts: outcome noise, estimated by nearest-neighbour
differencing in the estimator's own metric, and the variation from which units
are drawn, under the fixed-group-size design. Coverage of 95% intervals for
τ_super (`results/stage4_report.md`, `figures/coverage_vs_paper.png`):

| | Sim 1 ours | Sim 1 paper | Sim 2 ours | Sim 2 paper | SE / Monte Carlo SD (Sim 1, Sim 2) |
|---|---|---|---|---|---|
| CM | 94.5 [92.0, 97.2] | 95.0 [92, 98] | 90.3 [63.4, 98.0] | 97.1 [92, 100] | 1.05, 1.14 |
| GPS | 99.8 [99.6, 100] | 99.7 [99, 100] | 99.7 [98.6, 100] | 98.9 [96, 100] | 1.58, 1.51 |
| GPS_true | 94.4 | – | 94.2 | – | 0.98, 1.00 |
| DIF | 86.9 | – | 79.8 | – | 0.98, 1.01 |

Values are mean [min, max] over pairs, in %. The Monte Carlo SE of a coverage
near 95% is about 1 point.

* **GPS's over-coverage is reproduced, and explained.** The SE is calibrated
  when matching on the true scores (ratio 1.0, coverage 94%). With estimated
  scores it is about 1.5× too large, because the estimator then varies much
  less (SD 0.76–0.79 vs 1.23–1.25 with the true scores). The variance formula,
  like any that ignores score estimation, misses that reduction. Coverage of
  99.7–99.8% matches the paper's 99–100%.
* **CM matches the paper in Sim 1** (94.5% vs 95.0%). In Sim 2 it under-covers
  for the pairs with the largest bias: per-pair coverage correlates −0.92 with
  |bias|/SD, and even the oracle interval built from the true SD covers only
  85.6% on average. This follows from the Stage 3 finding that our CM is more
  biased than the paper's.
* DIF under-covers because of its confounding bias. Its SE is calibrated.

## Stage 1 findings (details in `docs/decisions.md`)

* $Z_1$ ("N[−3,3]") is read as Uniform(−3, 3). The centring in Sim 2 then makes
  every true effect zero (D1).
* Each treatment group is sampled to exactly 100 (Sim 1) or 200 (Sim 2) units
  from the assignment model. Bias is measured against the population effect (D3,
  D4).
* The paper's signed biases appear reversed relative to its own definition (D5).
* The paper's DIF "MSE" values cannot be MSEs of a difference in means (D15).
* Some differences in the paper's DIF column remain unexplained (D20).
