# Comparison interface: rebuild 2

This document describes rebuild 2's results in a fixed form, so that a later
comparison with other implementations can use them without reading the code.
It covers **only rebuild 2**. It was written before any comparison, without
access to any other rebuild or cross-check of this paper, and it makes no
assumptions about them.

## Identity

| Item | Value |
|---|---|
| Repository | `wangqiu-syr/glw2022-rebuild2` |
| Version for comparison | 1.0.0 (tag `v1.0.0` on `main`) |
| Environment | R 4.3.3 (2024-02-29), Ubuntu 24.04.5 LTS, R from the Ubuntu archive |
| Packages | base R for all estimators; testthat 3.2.1 (tests), nnet 7.3-19 (one reference test), MASS 7.3-60.0.1 (Stage-1 checks), ggplot2 3.4.4 (figures) |
| Replications | 500 per simulation |
| Group sizes | exactly 100 (Sim 1, n = 400) or 200 (Sim 2, n = 1800) per treatment combination |
| Random numbers | L'Ecuyer-CMRG; master seed 20220401 (Sim 1) and 20220402 (Sim 2); stream r = `parallel::nextRNGStream` applied r times; truths from 10^8 draws, seed 777 |

## Definitions

- **Pairs.** $\tau(t,t') = E[Y(t) - Y(t')]$ for $t < t'$, in the order (1,2), (1,3),
  …, i.e. the order of the paper's Tables 1 and 3.
- **Bias.** $E(\hat\tau) - \tau$, the mean over replications minus the truth. The
  paper's printed signs are the reverse (decisions D5).
- **Truths** (column `truth`):
  - `tau_super`: the population effect under the stated DGP, computed
    analytically. Sim 1: −1.5, −1, −1.5, 0.5, 0, −0.5. Sim 2: 0 for all pairs.
    This is the primary truth and the one used against the paper.
  - `tau_mix`: the effect over the equal-weight mixture of treatment groups that
    the fixed-group-size sample represents. Monte Carlo error ≤ 0.001.
  - `sample_effect`: the replication's own sample-average effect (used for
    conditional coverage only).
- **Intervals.** Estimate ± 1.96 × SE, with the analytic variance of D26. DIF uses
  $s_t^2/n_t + s_{t'}^2/n_{t'}$. The conditional interval uses the noise part
  only.

## Estimator labels

| Label | Meaning |
|---|---|
| `DIF` | difference of group means |
| `CM` | covariates matching: 1-NN with replacement, Mahalanobis(X) + Mahalanobis(Z), pooled covariances |
| `GPS` | GPS matching on the scalar $\hat p_t$ from a main-effects multinomial logit (primary) |
| `GPS_cs` | as GPS, with a correctly specified logit (adds $x_1^2..x_4^2, z_1^2, z_2^2$) |
| `GPS_true` | as GPS, on the true propensity scores |
| `CM_stacked` | CM with one Mahalanobis distance on (X, Z) |
| `CM_xmahal_zeuclid` | CM with Mahalanobis(X) + Euclidean(Z) |

`CM_stacked` and `CM_xmahal_zeuclid` have bias and MSE only, no intervals.

## The comparison file

`results/rebuild2_summary_long.csv`: one row per simulation × estimator ×
pair × measure × truth. It has 3,570 rows and is rebuilt by
`scripts/export_summary_long.R` from tracked tables.

| Column | Content |
|---|---|
| `sim` | `sim1` or `sim2` |
| `estimator` | a label from the table above |
| `t`, `tprime` | the pair, $t < t'$ |
| `measure` | see below |
| `truth` | `tau_super`, `tau_mix`, `sample_effect`, or empty if the measure does not use a truth |
| `tau` | the truth's value (empty for `sample_effect`, which varies by replication) |
| `value` | the estimate of the measure |
| `mcse` | its Monte Carlo SE (empty for `se_ratio`) |
| `reps` | replications (500) |

| `measure` | Definition | MC SE |
|---|---|---|
| `mean_est` | mean of $\hat\tau$ over replications | SD/√R |
| `bias` | `mean_est` − τ | SD/√R |
| `abs_bias` | \|bias\| | SD/√R |
| `mse` | mean of $(\hat\tau - \tau)^2$ | SD of squared errors/√R |
| `sd` | SD of $\hat\tau$ over replications | SD/√(2(R − 1)) |
| `coverage` | share of 95% intervals containing τ | √(p(1 − p)/R) |
| `coverage_cond` | share of conditional 95% intervals containing the replication's sample-average effect | √(p(1 − p)/R) |
| `coverage_oracle` | share of $\hat\tau$ ± 1.96 × `sd` containing `tau_super` (diagnostic, not a usable method) | √(p(1 − p)/R) |
| `mean_se` | mean estimated SE | SD of SE/√R |
| `se_ratio` | `mean_se` / `sd` | – |

Coverage values are proportions (0–1), not percentages. Per-replication values
are in `results/<sim>/estimates.csv`, `results/<sim>/intervals.csv.gz` and
`results/cm_metric_sensitivity/estimates_<sim>.csv.gz`.

## Interpretive choices (crosswalk)

The reasons are in [`decisions.md`](decisions.md).

| ID | Point | Rebuild 2's choice |
|---|---|---|
| D1 | $Z_1$ ("N[−3,3]") | Uniform(−3, 3) |
| D2 | probabilities don't sum to 1 | standard softmax, $\eta_1 = 0$ |
| D3 | group sizes vs model-based assignment | exactly $n_t$ per group, drawn from the assignment model |
| D4 | truth | `tau_super` (primary), `tau_mix` (secondary) |
| D5 | bias sign | $E(\hat\tau) - \tau$; paper's signs reversed |
| D6 | Sim 2 $Y(1)$ lacks ε | ε added |
| D7 | Table 2 numbering; no clusters | one $A$-level factor |
| D8 | errors across potential outcomes | one $N(0,1)$ draw per unit |
| D9 | CM metric | Mahalanobis(X) + Mahalanobis(Z), pooled, divisor $n$ |
| D10 | matching | 1-NN, with replacement, ties averaged, own observation for own level |
| D11 | averaging population | all $n$ units |
| D12 | propensity model | main-effects multinomial logit, refitted in each replication |
| D13 | GPS score | scalar $\hat p_t$ for target level $t$ |
| D14 | trimming | none |
| D15 | MSE | standard; the paper's DIF "MSE" flagged as not an MSE |
| D17 | replications, random numbers | 500; one stream per replication |
| D21 | Sim 1 outcomes uncentred | as written |
| D22 | logit fitting | Newton–Raphson, converged in every replication |
| D23 | draw order | X via Cholesky, Z₁, Z₂, Z₃, then T by inverse CDF; batches of 1000; ε last |
| D24 | ties | relative tolerance 1e-10 |
| D25 | GPS sensitivity | `GPS_cs`, `GPS_true` |
| D26 | intervals | analytic variance, normal intervals |

## Suggested agreement rule

Two implementations agree on a measure for a pair when

$$|a - b| \le 2\sqrt{\mathrm{SE}_a^2 + \mathrm{SE}_b^2},$$

using each side's Monte Carlo SE. Results from different implementations or R
versions will not match bit for bit (random numbers, tie handling), so
agreement should be judged at the level of Monte Carlo error. The paper reports
no Monte Carlo SEs, so comparisons with the paper are descriptive.
