# Specification: Simulations 1 and 2 of Guo, Liu & Wang (2022)

Source: Guo, S., Liu, J., & Wang, Q. (2022). Effective learning during COVID-19:
Multilevel covariates matching and propensity score matching. *Annals of Data
Science, 9*(5), 967–982. Page numbers below are the journal's (967–982).

This specification was written from the paper only. Where the paper is unclear
or inconsistent, the choice made is stated here and argued in
[`decisions.md`](decisions.md) (cited as **D1**, **D2**, …). Section 5 of the
paper (pp. 973–975) defines the simulations; Sections 2–4 (pp. 970–973) define
the estimators.

---

## 1. Notation

| Symbol | Meaning |
|---|---|
| $X = (X_1,\dots,X_4)$ | individual-level covariates |
| $Z = (Z_1, Z_2, Z_3)$ | "school-level" covariates (in the simulations these are drawn per unit, like $X$; there is no cluster structure, **D7**) |
| $W = (X, Z)$ | all 7 covariates |
| $T \in \{1,\dots,A\}$ | treatment *combination* (level-1 option × level-2 option), $A=4$ in Sim 1 and $A=9$ in Sim 2 |
| $Y(t)$ | potential outcome under combination $t$; observed $Y_i = Y_i(T_i)$ |
| $n_t$, $n$ | group size and total sample size |
| $p_t(w) = \Pr(T = t \mid W = w)$ | generalized propensity score (GPS) |

Estimand (eq. 1, p. 970): $\tau(t,t') = E[Y(t) - Y(t')]$ for every pair $t < t'$:
6 pairs in Sim 1 and 36 in Sim 2.

---

## 2. Data-generating process

### 2.1 Covariates (both simulations)

* $X \sim N_4(0, \Sigma_X)$ with variances $(2,1,1,2)$ and
  $\mathrm{cov}(X_1,X_2)=1.0$, $\mathrm{cov}(X_1,X_3)=-1.0$, $\mathrm{cov}(X_1,X_4)=0.5$,
  $\mathrm{cov}(X_2,X_3)=-0.5$, $\mathrm{cov}(X_2,X_4)=0.2$, $\mathrm{cov}(X_3,X_4)=-0.3$:

  $$\Sigma_X=\begin{pmatrix} 2 & 1 & -1 & 0.5\\ 1 & 1 & -0.5 & 0.2\\ -1 & -0.5 & 1 & -0.3\\ 0.5 & 0.2 & -0.3 & 2\end{pmatrix}$$

  (positive definite; eigenvalues 3.444, 1.747, 0.496, 0.312).
* $Z_1 \sim \mathrm{Uniform}(-3, 3)$. The paper writes "$Z_1 \sim N[-3,3]$"
  (**D1**).
* $Z_2 \sim \chi^2_1$.
* $Z_3 \sim \mathrm{Bernoulli}(0.5)$.
* $X$, $Z_1$, $Z_2$, $Z_3$ are mutually independent, and units are i.i.d.

So $E[W] = (0,0,0,0,\,0,\,1,\,0.5)$.

### 2.2 Treatment assignment

Multinomial logit with baseline $\eta_1 \equiv 0$:

$$p_t(w) = \frac{\exp\{\eta_t(w)\}}{\sum_{s=1}^{A}\exp\{\eta_s(w)\}}.$$

The paper's denominator, $1+\sum_{t=1}^{A}\exp\{\eta_t\}$, does not sum to one
(**D2**).

**Simulation 1** ($A = 4$, p. 974):

| $t$ | level-1, level-2 | $\eta_t(x,z)$ |
|---|---|---|
| 1 | 0, 0 | $0$ |
| 2 | 0, 1 | $0.1(x_1^2 + x_2 + x_3 + x_4 + z_1 + z_2 + z_3)$ |
| 3 | 1, 0 | $0.1(x_1 + x_2^2 + x_3^2 + x_4 + z_1 + z_2 + z_3)$ |
| 4 | 1, 1 | $0.1(x_1 + x_2^2 + x_3^2 + x_4^2 + z_1 + z_2 + z_3)$ |

**Simulation 2** ($A = 9$, p. 974). Here $\eta_1,\dots,\eta_4$ are the same as in Sim 1, and:

| $t$ | $\eta_t(x,z)$ |
|---|---|
| 1 | $0$ |
| 2 | $0.1(x_1^2 + x_2 + x_3 + x_4 + z_1 + z_2 + z_3)$ |
| 3 | $0.1(x_1 + x_2^2 + x_3^2 + x_4 + z_1 + z_2 + z_3)$ |
| 4 | $0.1(x_1 + x_2^2 + x_3^2 + x_4^2 + z_1 + z_2 + z_3)$ |
| 5 | $0.1(x_1 + x_2^2 + x_3^2 + x_4^2 + z_1^2 + z_2 + z_3)$ |
| 6 | $0.1(x_1 + x_2^2 + x_3^2 + x_4^2 + z_1^2 + z_2^2 + z_3)$ |
| 7 | $0.1(x_1 + x_2^2 + x_3^2 + x_4^2 + z_1^2 + z_2^2 + z_3^2)$ |
| 8 | $0.1(x_1^2 + x_2^2 + x_3^2 + x_4^2 + z_1^2 + z_2^2 + z_3^2)$ |
| 9 | $0.1(x_1 + x_2 + x_3 + x_4 + z_1 + z_2 + z_3)$ |

($z_3^2 = z_3$ because $Z_3$ is binary, so $\eta_6 \equiv \eta_7$. Treatments 6 and 7
therefore have identical propensity functions and identical covariate
distributions given $T$. They differ only in their outcome models.) Under 2.1, the marginal treatment shares
$\pi_t = E[p_t(W)]$ (computed by Monte Carlo integration, `results/stage1/`) are:
Sim 1: 0.183 / 0.255 / 0.252 / 0.310; Sim 2: 0.060 / 0.082 / 0.081 / 0.097 /
0.129 / 0.148 / 0.148 / 0.186 / 0.068. The $z_2^2$ terms in $\eta_6,\eta_7,\eta_8$
($Z_2^2 = N(0,1)^4$) give strong separation in the tail of $Z_2$, which the paper
mentions on p. 975.

Table 2 of the paper labels the nine combinations as level-1 group (A, B, C) ×
level-2 cohort (a, b, c), but the numbering it prints is inconsistent (**D7**).
The DGP is defined directly by the index $t = 1,\dots,9$, so the labels do not
affect anything.

### 2.3 Outcomes

$\varepsilon \sim N(0,1)$, independent of everything else, with one draw per unit
for the observed outcome (**D8**).

**Simulation 1** (p. 974): $Y(t) = [x, z]^\top \beta_t + \varepsilon$, with no intercept, where

| $t$ | $\beta_t$ on $(x_1,x_2,x_3,x_4,z_1,z_2,z_3)$ |
|---|---|
| 1 | (1, 1, 1, 1, 1, 1, 1) |
| 2 | (2, 3, 1, 2, 2, 2, 2) |
| 3 | (3, 1, 2, 1, 3, 1, 3) |
| 4 | (4, 1, 3, 4, 4, 1, 4) |

**Simulation 2** (p. 975): $Y(t) = \beta_t^\top (x_1,x_2,x_3,x_4,\,z_1,\,z_2-1,\,z_3-0.5) + \varepsilon$, where

| $t$ | coefficients on $(x_1,x_2,x_3,x_4,z_1,(z_2-1),(z_3-0.5))$ |
|---|---|
| 1 | (1, 1, 1, 1, 1, 1, 1) |
| 2 | (2, 3, 1, 2, 2, 2, 2) |
| 3 | (3, 1, 2, 1, 3, 1, 3) |
| 4 | (4, 1, 3, 4, 4, 1, 4) |
| 5 | (5, 1, 1, 1, 1, 1, 1) |
| 6 | (6, 3, 1, 2, 2, 2, 9) |
| 7 | (7, 1, 8, 1, 3, 1, 8) |
| 8 | (8, 1, 5, 2.5, 4, 1, 7) |
| 9 | (9, 1, 7, 1.5, 4, 1, 7) |

Rows 1–4 are Sim 1's $\beta_1,\dots,\beta_4$ with $Z_2$ and $Z_3$ centred. The paper
prints $Y(1)$ without $+\varepsilon$; we add it (**D6**).

### 2.4 Sample size and sampling design (**D3**)

* Sim 1: $n_t = 100$ for every $t$, so $n = 400$.
* Sim 2: $n_t = 200$ for every $t$, so $n = 1800$.

To get exactly $n_t$ units per group while treatment is still assigned by the
model in 2.2, each replication draws candidate units $(W, T)$ i.i.d. from 2.1–2.2
in fixed-size batches. A candidate is kept if its group still has fewer than $n_t$
members and discarded otherwise, and drawing stops when every group is full. This
is the same as drawing $n_t$ i.i.d. units from $W \mid T = t$ for each $t$,
independently across groups. Outcomes are then generated for the kept units.

---

## 3. Estimators

All three estimators take the form $\hat\tau(t,t') = \hat m_t - \hat m_{t'}$, where
$\hat m_t$ depends on $t$ only. So each estimator's bias satisfies
$\mathrm{Bias}(t,t') = h_t - h_{t'}$. This holds for the paper's Table 3 to within
rounding (Section 7).

### 3.1 DIF: the sample difference

$$\hat\tau_{DIF}(t,t') = \bar Y_t - \bar Y_{t'}, \qquad \bar Y_t = \frac{1}{n_t}\sum_{i:T_i=t} Y_i .$$

### 3.2 CM: covariates matching (Section 3, p. 971, eq. 2)

For each unit $i$ ($i = 1,\dots,n$, in all groups) and each level $t$:

$$j^*(i,t) = \arg\min_{j:T_j=t}\big\{\|X_j - X_i\|_{V_X} + \|Z_j - Z_i\|_{V_Z}\big\}, \qquad
\|a\|_{V} = (a^\top V^{-1} a)^{1/2},$$

where $V_X = \frac1n\sum_{i=1}^n (X_i-\bar X)(X_i-\bar X)^\top$ over the whole sample,
and $V_Z$ is defined the same way for $Z$ (**D9**). The imputed potential outcome is
$\hat Y_i(t) = Y_{j^*(i,t)}$. For $T_i = t$ this is unit $i$ itself (distance 0),
so $\hat Y_i(t) = Y_i$. Matching uses one nearest neighbour, with replacement; ties
are averaged (**D10**). Then

$$\hat\tau_{CM}(t,t') = \frac1n\sum_{i=1}^n\big\{\hat Y_i(t) - \hat Y_i(t')\big\},$$

averaged over all $n$ units, including those in neither group $t$ nor group $t'$ (**D11**).

### 3.3 GPS: generalized propensity score matching (Section 4, p. 972, eqs. 3–4)

1. In each replication, fit a multinomial logistic regression of $T$ on an
   intercept and the main effects $X_1,\dots,X_4,Z_1,Z_2,Z_3$ by maximum
   likelihood. This model is misspecified, because the true $\eta_t$ contain
   squares (**D12**). It gives $\hat p_t(W_i)$ for all $i$ and $t$.
2. For each unit $i$ and level $t$, match on the scalar score for level $t$:
   $j^*(i,t) = \arg\min_{j:T_j=t} |\hat p_t(W_j) - \hat p_t(W_i)|$, with one
   neighbour, replacement, and own observation when $T_i = t$ (**D10**, **D13**).
3. $\hat Y_i(t) = Y_{j^*(i,t)}$ and
   $\hat\tau_{GPS}(t,t') = \frac1n\sum_{i=1}^n\{\hat Y_i(t) - \hat Y_i(t')\}$.

No trimming is applied (**D14**).

Sensitivity variants (**D25**), computed in every replication with the same
matching step: **GPS_cs** uses a correctly specified multinomial logit (main
effects plus $x_1^2..x_4^2, z_1^2, z_2^2$), and **GPS_true** uses the true $p_t(W)$.

### 3.4 Not part of the reproduction

The covariate-balance algorithms (Steps 1–3 on pp. 971 and 973) are not used in
any simulation table. They may be implemented later as diagnostics; their
inconsistencies are listed in **D18**. Section 6 (synthetic data analysis) is
outside scope.

---

## 4. True effects

**Primary truth: the population effect $\tau_{super}$** (eq. 1 under the DGP of
Section 2, **D4**). The outcome models are linear in $W$, so
$\tau(t,t') = E[W]^\top(\beta_t - \beta_{t'})$ plus the difference in intercepts:

* **Sim 1**: $E[Y(t)] = 1.5,\ 3,\ 2.5,\ 3$ for $t = 1,\dots,4$, which gives

  | pair | (1,2) | (1,3) | (1,4) | (2,3) | (2,4) | (3,4) |
  |---|---|---|---|---|---|---|
  | $\tau$ | −1.5 | −1.0 | −1.5 | 0.5 | 0.0 | −0.5 |

* **Sim 2**: the centring makes $E[Y(t)] = 0$ for every $t$, so $\tau(t,t') = 0$
  for all 36 pairs.

**Secondary truth: the mixture effect $\tau_{mix}$.** Under design 2.4 the pooled
sample is an equal-weight mixture of the group-conditional distributions, so
CM and GPS, which average over all $n$ sampled units, target

$$\tau_{mix}(t,t') = \frac1A\sum_{s=1}^{A} E\big[\mu_t(W) - \mu_{t'}(W) \mid T = s\big], \qquad \mu_t(w) = E[Y(t)\mid W=w].$$

It is computed by high-precision Monte Carlo integration. Stage 1 values (N = 4×10⁶)
differ from $\tau_{super}$ by at most about 0.09 (Sim 1) and 0.17 (Sim 2). Bias
is reported against both truths, and against $\tau_{super}$ when compared with the
paper.

---

## 5. Performance measures

For each simulation, estimator $e \in \{\mathrm{DIF}, \mathrm{CM}, \mathrm{GPS}\}$ and
pair $(t,t')$, over $R$ replications with estimates $\hat\tau_r$:

| Measure | Definition | In paper |
|---|---|---|
| Bias | $\bar{\hat\tau} - \tau$, with $\bar{\hat\tau} = R^{-1}\sum_r \hat\tau_r$ (sign convention: **D5**) | Table 3 (signed) |
| Absolute bias | $\lvert\bar{\hat\tau} - \tau\rvert$ (Table 1 footnote: "$\lvert E(\hat\tau) - \tau\rvert$") | Table 1 |
| MSE | $R^{-1}\sum_r(\hat\tau_r - \tau)^2$ (**D15**) | Table 3 |
| Empirical SD | $\{(R-1)^{-1}\sum_r(\hat\tau_r - \bar{\hat\tau})^2\}^{1/2}$ | no (added) |
| Monte Carlo SEs | of bias ($\mathrm{SD}/\sqrt R$) and of MSE | no (added) |
| Coverage | share of replications whose 95% interval contains $\tau$; interval method: analytic variance with normal intervals (**D26**) | Tables 1 and 3, CM and GPS only, whole percent |

What the paper reports:

* **Table 1** (Sim 1, "500 Monte Carlo datasets"): absolute bias for DIF, CM and GPS;
  coverage (%) for CM and GPS.
* **Table 3** (Sim 2): bias and MSE for DIF, CM and GPS; coverage (%) for CM and GPS.
* **Fig. 1**: Table 3's bias columns plotted against pair index 1–36. It adds no new numbers.

We report bias, absolute bias and MSE (plus SD and MC SEs) for all estimators
and pairs in both simulations. Sim 1's signed bias and MSE have no counterpart in
the paper.

---

## 6. Replications and random numbers

* $R = 500$ replications for each simulation (p. 974: "over 500 Monte Carlo datasets";
  p. 975: "We repeat the experiment 500 times").
* One independent random-number stream per replication: R's `L'Ecuyer-CMRG`
  generator, where stream $r$ is obtained by applying `parallel::nextRNGStream`
  $r$ times to a fixed master seed. Replication $r$ is therefore reproducible
  without running the others. Each replication's estimates are saved to disk as
  soon as it finishes, so a run can resume.
* Software: R 4.3.3 (Ubuntu 24.04 archive). The code in `R/` uses base R only
  (including `parallel`). The multinomial logit is a hand-written Newton–Raphson
  (**D22**). Packages: testthat 3.2.1 (tests); nnet 7.3-19 (only as a reference
  in one test of the logit fit); MASS 7.3-60.0.1 (only in the Stage-1 diagnostic
  scripts); ggplot2 3.4.4 (figures, Stage 3). No package implementing a matching
  estimator is used.

---

## 7. The paper's tables as CSV, and checks of the transcription

Files (values exactly as printed; coverage in percent, without the "%"):

* `data/paper/table1_sim1.csv`: `effect, t, tprime, dif_absbias, cm_absbias, cm_coverage, gps_absbias, gps_coverage` (6 rows)
* `data/paper/table3_sim2.csv`: `effect, t, tprime, dif_bias, dif_mse, cm_bias, cm_mse, cm_coverage, gps_bias, gps_mse, gps_coverage` (36 rows)

Checks:

1. **Text layer compared with the rendered page.** The values were parsed from
   the PDF text layer (`pdftotext -layout`), then compared cell by cell with
   300-dpi renders of pp. 974, 976 and 977. All 6 + 36 rows agree in every column.
2. **Internal consistency** (`scripts/stage1_check_tables.R`). Each bias column
   must have the form $h_t - h_{t'}$ (Section 3). Fitting 8 free parameters
   to 36 pairs gives a maximum absolute residual of 0.00006 (DIF), 0.00007 (CM)
   and 0.00006 (GPS), i.e. rounding error at 4 decimals. For Table 1, which
   gives absolute values only, a sign pattern exists that makes each column
   exactly additive (residual 0.0000): $(+,+,+,-,+,+)$ for DIF over pairs
   (1,2), (1,3), (1,4), (2,3), (2,4), (3,4), and all signs equal for CM and for GPS.
   This also confirms that the paper's estimators impute $\hat Y_i(t)$
   independently of $t'$, as specified in Section 3.

The per-treatment decomposition $h_t$ of the paper's Table 3 biases is saved in
`results/stage1/paper_table3_bias_decomposition.csv`.

---

## 8. Stage-1 diagnostic evidence used for decisions

The DIF estimator involves no modelling choices, so its expected value depends
only on the DGP. Its population bias can be computed by Monte Carlo integration
(`scripts/stage1_dgp_checks.R`, weighting by the true $p_t$ with N = 4×10⁶) and
compared with the paper's DIF column. The paper's Monte Carlo error in that
column is about 0.02–0.06. Summary (details in `results/stage1/` and **D1–D5**):

| Check | Result |
|---|---|
| Sim 2, $Z_1$ variants: RMS difference from the paper's 36 DIF biases (sign-aligned) | U(−3,3): 0.40 (cor 0.956); N(0, var 3): 0.46; N(−3, var 3): 2.36; N(−3, sd 3): 5.07 |
| Sim 2, sign of bias | Correlation with the paper is −0.956 as defined, +0.956 with the sign flipped (**D5**) |
| Normalisation variants (Sim 2 RMS) | softmax 0.403; leftover mass to $t=1$: 0.416; to $t=A$: 0.792 |
| Paper's CM bias pattern (Sim 2) | Correlation 0.82 (30/36 signs agree) with the sign-aligned $\tau_{mix}-\tau_{super}$, against 0.44 with the population DIF bias. This supports design 2.4 together with truth $\tau_{super}$ (**D3**, **D4**) |
| Paper's DIF "MSE" (Sim 2) | 49–335, about 0.81–1.00 × $\{\mathrm{Var}(Y\mid T=t)+\mathrm{Var}(Y\mid T=t')\}$; closer to these variances than to an MSE of a difference in means (**D15**) |
| Remaining DIF mismatch | Sim 1: pairs involving $t=3$ (population $h_3$ ≈ 0.70 away from the paper's); Sim 2: $t=4$ and $t=7$ (≈ 0.4–0.6). Not explained by any variant tried; flagged in **D20** |
