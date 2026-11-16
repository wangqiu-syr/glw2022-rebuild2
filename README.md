# glw2022-rebuild2

An independent rebuild, in R, of the two simulation studies in

> Guo, S., Liu, J., & Wang, Q. (2022). Effective learning during COVID-19:
> Multilevel covariates matching and propensity score matching.
> *Annals of Data Science, 9*(5), 967–982.

It was built from the paper alone. No code, replications or discussion of the
paper were used, nor the cited works' implementations. The matching estimators
are written from scratch. The paper's PDF is not in this repository (copyright).

## Status

| Stage | Content | Status |
|---|---|---|
| 1 | Specification, decisions, paper's tables as CSV | **done, awaiting review** |
| 2 | Code in `R/`, tests, pilot run | not started |
| 3 | Full runs (500 replications each): bias, \|bias\|, MSE against the paper | not started |
| 4 | Intervals and coverage (optional) | not started |

## Stage 1 contents

* [`docs/spec.md`](docs/spec.md): the DGP, the estimators (DIF, CM, GPS), the true
  effects, the performance measures, the replications, and the stage-1 evidence.
* [`docs/decisions.md`](docs/decisions.md): 21 points where the paper is unclear or
  inconsistent, each with options, choice and reason.
* [`data/paper/`](data/paper): Tables 1 and 3 as CSV, checked against the
  rendered PDF and for internal consistency.
* `scripts/stage1_*.R` and `results/stage1/`: population-level checks used to
  settle the unclear DGP details, with their outputs.

Main stage-1 findings:
* $Z_1$ ("N[−3,3]") is read as Uniform(−3, 3). Among the readings tried, the
  population DIF bias then fits the paper's DIF column best, and the centring
  in Sim 2 makes every true effect zero (D1).
* Each treatment group is sampled to exactly 100 (Sim 1) or 200 (Sim 2) units
  from the assignment model. Bias is measured against the population effect. The
  paper's CM biases carry the signature of this design (D3, D4).
* The paper's signed biases appear reversed relative to its own definition (D5).
* The paper's DIF "MSE" values (49–335) cannot be MSEs of a difference in means.
  They are close to the sum of the two within-group outcome variances (D15).
* Some differences in the DIF column remain unexplained (D20).

## Reproducing stage 1

R 4.3.3 with MASS, from the Ubuntu 24.04 archive (`r-base-core`, `r-cran-mass`).
From the repository root:

```sh
Rscript scripts/stage1_check_tables.R     # transcription consistency (seconds)
Rscript scripts/stage1_dgp_checks.R       # population checks, Z1 variants (~2.5 min)
Rscript scripts/stage1_dgp_variants.R     # normalisation variants (~1 min)
Rscript scripts/stage1_design_evidence.R  # design/truth evidence (seconds)
```
