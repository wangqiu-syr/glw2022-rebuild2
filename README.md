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
| 2 | Code in `R/`, tests, population truths, pilot run | **done, awaiting review** |
| 3 | Full runs (500 replications each): bias, \|bias\|, MSE against the paper | not started |
| 4 | Intervals and coverage (optional) | not started |

## Layout

| Path | Content |
|---|---|
| [`docs/spec.md`](docs/spec.md) | DGP, estimators, true effects, performance measures, replications |
| [`docs/decisions.md`](docs/decisions.md) | 25 recorded choices where the paper is unclear or inconsistent, or where implementation needed a choice |
| [`data/paper/`](data/paper) | the paper's Tables 1 and 3 as CSV (checked) |
| `R/config.R` | simulation settings (A, n per group, master seeds, estimators) |
| `R/dgp.R` | covariates, treatment model, outcomes, group-filling sampler |
| `R/mlogit.R` | multinomial logit by Newton–Raphson; GPS design matrices |
| `R/estimators.R` | DIF, CM (Mahalanobis 1-NN), GPS matching |
| `R/truth.R` | τ_super (analytic), τ_mix and other population quantities (Monte Carlo integration) |
| `R/rng.R`, `R/run.R` | per-replication random streams; resumable parallel runner; result collection |
| `R/performance.R` | bias, \|bias\|, MSE, SD, Monte Carlo SEs |
| `tests/testthat/` | 87 tests, including hand-worked examples |
| `scripts/` | drivers: tests, truths, simulation runs, pilot summary, stage-1 checks |
| `results/truth/` | population truths per simulation |
| `results/pilot/` | 20-replication pilot: per-replication files, performance preview, timing |
| `results/stage1/` | stage-1 diagnostic outputs |

Estimators computed in every replication: `DIF`, `CM`, `GPS` (main-effects
logit, primary), plus the sensitivity variants `GPS_cs` (correctly specified
logit) and `GPS_true` (true propensity scores).

## Reproducing

R 4.3.3 from the Ubuntu 24.04 archive (`r-base-core`, `r-cran-testthat`; `r-cran-nnet`
and `r-cran-mass` for one test and the stage-1 scripts). From the repository root:

```sh
Rscript scripts/run_tests.R                                   # 87 tests (~5 s)
Rscript scripts/compute_truth.R                               # population truths (~2.5 min, 4 cores)
Rscript scripts/run_simulations.R sim1,sim2 1 20 results/pilot/runs 4   # pilot (~10 s)
Rscript scripts/pilot_summary.R                               # results/pilot/pilot_summary.md
```

`run_simulations.R` saves each replication as it finishes and skips finished ones
when rerun, so an interrupted run can be resumed. Replication *r* uses its own
L'Ecuyer-CMRG stream, so it reproduces whether run alone or in a batch.

Stage-1 checks (population-level comparison with the paper's DIF column):

```sh
Rscript scripts/stage1_check_tables.R; Rscript scripts/stage1_dgp_checks.R
Rscript scripts/stage1_dgp_variants.R; Rscript scripts/stage1_design_evidence.R
```

## Findings so far

Stage 1 (details in `docs/decisions.md`):
* $Z_1$ ("N[−3,3]") is read as Uniform(−3, 3). The centring in Sim 2 then makes
  every true effect zero (D1).
* Each treatment group is sampled to exactly 100 (Sim 1) or 200 (Sim 2) units
  from the assignment model. Bias is measured against the population effect (D3,
  D4).
* The paper's signed biases appear reversed relative to its own definition (D5).
* The paper's DIF "MSE" values cannot be MSEs of a difference in means (D15).
* Some differences in the paper's DIF column remain unexplained (D20).

Stage 2 pilot (20 replications; `results/pilot/pilot_summary.md`): all logit fits
converged; the simulated DIF bias agrees with the population calculation;
projected wall-clock for 500 replications on 4 cores is about 0.3 min (Sim 1)
and 2.4 min (Sim 2).
