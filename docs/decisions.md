# Decisions log

For every point where Guo, Liu & Wang (2022) is unclear, ambiguous or
inconsistent: the issue, the options, the choice made, and the reason. Page
numbers are the journal's. Impact: **H** means it can change the reported
numbers materially, **M** means it changes them somewhat, **L** means cosmetic or
no effect.

The "Evidence" lines refer to Stage-1 checks that compare population quantities,
computed by Monte Carlo integration, with the paper's own DIF column. DIF is the
one estimator whose expectation depends only on the DGP (`scripts/stage1_*.R`,
`results/stage1/`).

Summary:

| ID | Topic | Choice | Impact |
|---|---|---|---|
| D1 | Distribution of $Z_1$ ("N[−3,3]") | Uniform(−3, 3) | H |
| D2 | Softmax denominator does not sum to 1 | Standard softmax, $\eta_1 = 0$ | L |
| D3 | "Sample size for each treatment combination is 100/200" vs model-based assignment | Exactly $n_t$ per group, drawn from the assignment model | H |
| D4 | Which "true" effect | Population $\tau_{super}$ (primary); $\tau_{mix}$ also reported | H |
| D5 | Sign convention of bias | $E(\hat\tau)-\tau$, $\tau(t,t')=E[Y(t)-Y(t')]$; paper's signed values compared after a sign flip | M |
| D6 | $Y(1)$ in Sim 2 has no $\varepsilon$ | Add $\varepsilon$ | L |
| D7 | Table 2 numbering; no cluster structure in the DGP | Treat $t$ as a 9-level factor; labels cosmetic | L |
| D8 | Error terms across potential outcomes | One $N(0,1)$ draw per unit | L |
| D9 | CM metric | Sum of separate Mahalanobis distances for $X$ and $Z$, full-sample covariances | M |
| D10 | Number of matches, replacement, ties, own observation | 1-NN, with replacement, ties averaged, own observation kept | M |
| D11 | Population over which CM/GPS average | All $n$ units | M |
| D12 | GPS model | Main-effects multinomial logit, fitted by ML in each replication | H |
| D13 | GPS matching score | Scalar $\hat p_t$ for target level $t$ | M |
| D14 | Trimming | None | M |
| D15 | Definition of MSE; paper's DIF MSE | Standard MSE; paper's DIF MSE flagged as a different quantity | H (for comparison) |
| D16 | Interval method for coverage | Deferred to Stage 4 | — |
| D17 | Replications and random numbers | 500 each; one L'Ecuyer-CMRG stream per replication | L |
| D18 | Balance-diagnostic formulas | Out of scope; typos noted | — |
| D19 | Extra measures not in the paper | Signed bias and MSE for Sim 1; empirical SD; MC SEs | — |
| D20 | Residual mismatch in the DIF column | Unresolved; flagged | — |
| D21 | Sim 1 outcomes not centred (Sim 2 are) | Taken literally | L |
| D22 | Fitting the multinomial logit (Stage 2) | Hand-written Newton–Raphson, checked against `glm` and `nnet::multinom` | L |
| D23 | Order and method of random draws (Stage 2) | X via Cholesky factor; T by inverse CDF; batches of 1000 | L |
| D24 | Tie tolerance in matching (Stage 2) | Relative 1e-10 | L |
| D25 | Sensitivity GPS models (agreed at Stage 1 review) | GPS_cs (correctly specified logit) and GPS_true (true scores) | — |

---

## D1. Distribution of $Z_1$ — impact H

**Issue.** p. 973: "$Z_1 \sim N[-3, 3]$". "N" usually means normal (with mean −3
and variance or SD 3?). The square brackets are interval notation, which suggests
a uniform on $[-3, 3]$.

**Options.** (a) Uniform(−3, 3); (b) $N(-3, 3)$, variance 3; (c) $N(-3, 3^2)$;
(d) $N(0, 3)$.

**Choice.** (a) Uniform(−3, 3).

**Reason.**
1. The bracketed interval reads naturally as a support.
2. In Sim 2 the authors centre $Z_2$ and $Z_3$ ($z_2 - 1$, $z_3 - 0.5$) but not
   $Z_1$. That only makes all $E[Y(t)]$ equal (to zero) if $E[Z_1] = 0$, which
   rules out the mean −3 readings.
3. *Evidence.* Root-mean-square difference between the population DIF bias and
   the paper's 36 DIF biases (sign-aligned, D5): U(−3,3) 0.40, N(0, var 3) 0.46,
   N(−3, var 3) 2.36, N(−3, sd 3) 5.07 (`results/stage1/dgp_population_checks_summary.csv`).
   The two mean-zero readings fit far better. Between them the uniform fits
   slightly better and matches the notation. Sim 1 hardly discriminates (RMS
   0.31, 0.31, 0.33, 0.84).

---

## D2. Normalisation of the multinomial-logit probabilities — impact L

**Issue.** pp. 973–974: $p_t = \exp\{\eta_t\}/[1+\sum_{t=1}^{A}\exp\{\eta_t\}]$
with $\eta_1 = 0$. The denominator then counts $\exp(\eta_1) = 1$ twice, so
$\sum_t p_t = 1 - 1/D < 1$. That contradicts the stated $\sum_t p_t = 1$.

**Options.** (a) Standard softmax, denominator $\sum_{s=1}^{A}\exp(\eta_s) = 1+\sum_{s\ge2}\exp(\eta_s)$;
(b) the literal formula, with the leftover mass assigned to $t = 1$; (c) leftover
assigned to $t = A$; (d) leftover units discarded, which is the same as (a).

**Choice.** (a).

**Reason.** It is the standard parameterisation, it satisfies the paper's stated
constraint, and any sampler that normalises (e.g. R's `sample(prob=)`) gives (a)
anyway. *Evidence* (`results/stage1/dgp_assignment_variants.csv`), RMS
against the paper's DIF bias: Sim 2: (a) 0.403, (b) 0.416, (c) 0.792; Sim 1:
(a) 0.309, (b) 0.360, (c) 0.594. (a) fits best in both.

---

## D3. Sampling design: fixed group sizes vs model-based assignment — impact H

**Issue.** Treatment is assigned by the multinomial model, which makes group
sizes random and unequal ($\pi_t$ ranges from 0.18 to 0.31 in Sim 1 and from 0.06
to 0.19 in Sim 2). But the paper says "the sample size for each treatment
combination is 100" (Sim 1, p. 974) and "the sample size of 1800 where each
treatment combination is 200 observations" (Sim 2, p. 975).

**Options.**
(a) Fixed total $n$ ($400$ or $1800$) with multinomial assignment, so group sizes
are random and equal 100/200 only on average, which they do not.
(b) Exactly $n_t$ per group. Draw candidate $(W, T)$ from the model and keep each
group's first $n_t$ members. This is equivalent to $n_t$ i.i.d. draws from
$W \mid T=t$.
(c) Draw $W$ for $n$ units and force 100/200 per group by another rule (e.g.
ranking on $p_t$). The paper does not describe any such rule.

**Choice.** (b) (spec §2.4).

**Reason.** It is the only option that satisfies both explicit statements: the
model-based assignment and equal group sizes, with $9 \times 200 = 1800$.
*Evidence.* Under (b) the pooled sample is an equal-weight mixture of the groups,
so CM and GPS target $\tau_{mix}$ rather than $\tau_{super}$ (D4). If the paper
measured bias against $\tau_{super}$, its CM biases should therefore track
$\tau_{mix} - \tau_{super}$. They do: correlation 0.82, with 30 of 36 signs
agreeing (sign-aligned). The correlation with the population DIF bias is only
0.44, which is what one would expect if CM bias were just leftover confounding
(`results/stage1/design_evidence_sim2.csv`). Option (a) gives no such shift.
Under (a) and (b) the DIF expectation is the same, so DIF cannot tell them apart.

---

## D4. The true effect $\tau$ — impact H

**Issue.** Eq. (1) defines $\tau(t,t') = E[Y(t)-Y(t')]$ but does not say over
which population, or how the simulation's truth was computed. Under D3(b) the
sample is not a random sample from the covariate distribution in spec §2.1.

**Options.** (a) $\tau_{super}$: the expectation under the covariate distribution
of §2.1, computed analytically. (b) $\tau_{mix}$: the expectation over the
equal-weight mixture of the group-conditional distributions, which is what the
sample represents under D3(b). (c) The per-replication sample average effect
$n^{-1}\sum_i\{\mu_t(W_i)-\mu_{t'}(W_i)\}$, whose expectation under D3(b) is $\tau_{mix}$.

**Choice.** (a) is primary and is used in all comparisons with the paper. (b) is
also reported, as the estimand CM and GPS actually target under D3(b).

**Reason.** Sim 2's careful centring of $Z_2$ and $Z_3$ only makes sense if the
truth comes from population means, giving $\tau = 0$ for all 36 pairs. The CM
evidence under D3 points the same way. Reporting (b) separates the estimators'
own bias from the shift in estimand caused by the design.
Values: Sim 1 $\tau_{super}$ = −1.5, −1, −1.5, 0.5, 0, −0.5 for pairs (1,2), (1,3),
(1,4), (2,3), (2,4), (3,4). Sim 2 $\tau_{super} = 0$. $\tau_{mix}$ will be computed
to high precision in Stage 2.

---

## D5. Sign convention of the bias — impact M (signed comparisons only)

**Issue.** Eq. (1) and the Table 1 footnote ("abs bias is $|E(\hat\tau)-\tau|$")
imply Bias $= E[\hat\tau(t,t')] - \tau(t,t')$ with $\tau(t,t') = E[Y(t) - Y(t')]$.
But the paper's signed DIF biases in Table 3 correlate −0.956 with the population
DIF bias computed under that convention, and +0.956 after flipping the sign.
So the paper appears to have used $\hat\tau(t,t')$ for $E[Y(t')-Y(t)]$, or
Bias $= \tau - E(\hat\tau)$. The two readings cannot be told apart.

**Options.** (a) Follow the written definitions and flip the paper's signed values
when comparing. (b) Adopt the paper's apparent convention.

**Choice.** (a). Our tables use Bias $= E(\hat\tau)-\tau$ with $\tau(t,t') = E[Y(t)-Y(t')]$.
In comparisons, the paper's signed Table 3 biases are multiplied by −1 and
labelled "paper, sign-aligned". Absolute bias and MSE do not depend on the sign.
For Table 1, which gives |bias| only, the internally consistent sign pattern (spec
§7) is used, with the same flip, wherever signed values are plotted.

**Reason.** The written definitions are explicit; the flip is an inference. It is
stated openly and is easy to reverse.

---

## D6. Missing error term in Sim 2's $Y(1)$ — impact L

**Issue.** p. 975: "$Y(1) = x_1 + x_2 + x_3 + x_4 + z_1 + z_2 - 1 + z_3 - 0.5$", with
no $+\epsilon$. Every other $Y(t)$ has one, as does Sim 1's $Y(t)$.

**Choice.** Add $\varepsilon \sim N(0,1)$.

**Reason.** Almost certainly a typo. Without it, only the variance of estimates
involving $t = 1$ would change slightly.

---

## D7. Table 2 numbering; no cluster structure — impact L

**Issue.** (i) Table 2 (p. 974) numbers the cells 1 = Aa, 2 = Ba, 3 = Ca,
**3** = Ab, **4** = Bb, 6 = Cb, 7 = Ac, 8 = Bc, 9 = Cc, so 3 appears twice and 5 is
missing. The intended numbering is presumably row-wise: Ab = 4, Bb = 5. (ii) The
"multilevel" structure (level-1 × level-2, school vs individual covariates)
has no role in the simulation DGP. $Z$ is drawn per unit, there are no clusters,
and both $\eta_t$ and $Y(t)$ are indexed by the combined $t$ directly.

**Choice.** Treat $T$ as a single factor with $A$ levels, as the DGP and the
estimators define it. Read Table 2 row-wise (1 = Aa, …, 9 = Cc) for labels only.

**Reason.** Nothing numerical depends on the labels.

---

## D8. Error terms across potential outcomes — impact L

**Issue.** The paper writes one $\epsilon$ in each $Y(t)$ without saying whether
the errors are shared or independent across $t$.

**Choice.** Draw one $\varepsilon_i \sim N(0,1)$ per unit and add it to the observed
$\mu_{T_i}(W_i)$.

**Reason.** Only the observed outcome is ever used, and the truth $\tau_{super}$
depends only on $\mu_t$. The joint distribution of unobserved errors therefore
affects nothing that is reported.

---

## D9. Metric for covariates matching — impact M

**Issue.** p. 971 defines $CM(t,x,z) = \arg\min_{j:T_j=t}(\|X_j-x\| + \|Z_j - z\|)$,
with $\|\cdot\|$ the Mahalanobis norm. The covariance $V = \frac1n\sum_i(X_i-\bar X)(X_i-\bar X)^\top$
is given for $X$ only.

**Options.** (a) Sum of two Mahalanobis distances, $\|\cdot\|_{V_X}$ on $X$ and
$\|\cdot\|_{V_Z}$ on $Z$, with $V_Z$ defined the same way. (b) One Mahalanobis
distance on the stacked $(X,Z)$. (c) Mahalanobis on $X$ plus Euclidean on $Z$.

**Choice.** (a). $V_X$ and $V_Z$ use the whole pooled sample with divisor $n$, as
written. $Z_3$ (binary) is included in the $Z$ block.

**Reason.** (a) follows the formula as written, two norms added, and $V_Z$
by analogy. The divisor ($n$ vs $n-1$) rescales both blocks by the same factor and
cannot change the argmin.

**Sensitivity check (after Stage 3, at your request).** CM was rerun under (b)
and (c) on the same 500 + 500 data sets (`results/cm_metric_sensitivity/`,
`figures/cm_metric_sensitivity_sim2.png`). The metric makes little difference.
Sim 2 mean |bias| is 0.299 (a), 0.304 (b) and 0.266 (c), against the paper's
0.139. The paper's CM bias is about a third of ours under every metric (slope
0.32–0.35), and mean MSE is 0.29 in all three. Sim 1 mean |bias| is 0.166, 0.179
and 0.146, against the paper's 0.143. So the metric does not explain why the
paper's CM is less biased. (a) stays the primary choice.

---

## D10. Matching details — impact M

**Issue.** The paper uses "argmin", which implies a single nearest neighbour.
It does not mention replacement, ties, or what happens when $T_i = t$.

**Choice.** One nearest neighbour ($M = 1$), with replacement (every unit in
group $t$ is a candidate for every $i$). When $T_i = t$ the argmin includes $j = i$
at distance 0, so $\hat Y_i(t) = Y_i$. Exact ties at the minimum distance, which
occur with probability 0 for continuous covariates except $j = i$, are resolved
by averaging the tied outcomes.

**Reason.** This is the literal reading of the argmin over $\{j : T_j = t\}$.

---

## D11. Population averaged over in CM and GPS — impact M

**Issue.** Eq. (2) and eq. (4) average over $i = 1,\dots,n$, and the paper notes
"We impute potential outcomes $Y_i$ for individuals who did not receive either
treatment $t$ or $t'$". So the average is over the full sample, not over the units
in groups $t$ and $t'$ only.

**Choice.** Average over all $n$ units, as written.

**Reason.** Explicit in the paper. Combined with D3(b), this is why $\tau_{mix}$
is the estimand that CM and GPS target (D4).

---

## D12. Specification of the GPS model — impact H

**Issue.** The paper never says how $\Pr(t\mid x,z)$ is estimated for GPS matching.
The true model has squared terms. The conclusion (p. 978) attributes GPS's lower
efficiency to "the misspecification of the propensity score modeling".

**Options.** (a) Multinomial logit with intercepts and main effects $X_1..X_4, Z_1..Z_3$,
which is misspecified. (b) The correctly specified model: linear plus the squared
terms that appear in some $\eta_t$, i.e. $x_1^2..x_4^2, z_1^2, z_2^2$. (c) The true
$p_t$.

**Choice.** (a), fitted by maximum likelihood in every replication. Possible
sensitivity analyses: (b) and (c), run only if you approve them at Stage 3.

**Reason.** The paper's conclusion implies a misspecified model, and (a) is the
default an analyst would fit without knowing the DGP. Software: `nnet::multinom`
(general model fitting, allowed by the brief) or a hand-written Newton–Raphson.
This is to be decided and tested in Stage 2, with convergence checked in every
replication.

---

## D13. Score used in GPS matching — impact M

**Issue.** Eq. (3) defines $GPS(t, p) = \arg\min_{j:T_j=t}\|\Pr(t\mid x_j,z_j) - p\|$ and
imputes $\hat Y_i(t) = Y^{obs}_{GPS}(t, \Pr(t\mid x_i,z_i))$.

**Choice.** For target level $t$, match unit $i$ to the unit in group $t$ whose
estimated $\hat p_t$ is closest to $\hat p_t(W_i)$, by absolute difference. A
different scalar is used for each target $t$. The full GPS vector is not used.

**Reason.** This is the literal reading, and the paper states that "the treatment
level $t$ enters the matching scheme through the function of the covariates that
is being matched, $\Pr(t\mid x_j,z_j)$".

---

## D14. Trimming — impact M

**Issue.** Section 2 (p. 970) says the authors "generalize the [47] approach [trimming]
to the multilevel treatment combinations via constructing subsamples", but gives
no rule, and Section 5 mentions no trimming.

**Choice.** No trimming. The full sample is used.

**Reason.** No rule is given that could be implemented. The paper's Table 3 bias
columns are exactly additive (spec §7), which is consistent with every unit
receiving every $\hat Y_i(t)$, i.e. with a fixed, untrimmed sample.

---

## D15. MSE definition, and the paper's DIF "MSE" — impact H for comparison

**Issue.** The paper's Table 3 DIF MSEs are 49–335. For a difference of two
group means with 200 units each, the MSE is
$\mathrm{bias}^2 + \{\mathrm{Var}(Y\mid T=t)+\mathrm{Var}(Y\mid T=t')\}/200$,
which is about 0.3–1.9 plus bias² (at most about 8 here). The paper's values are
instead 0.81–1.00 times $\mathrm{Var}(Y\mid T=t)+\mathrm{Var}(Y\mid T=t')$ itself,
computed under our DGP (`scripts/stage1_dgp_checks.R`, column `varsum`). That
would fit, for example, variances not divided by group size. The CM and GPS MSEs
(0.13–3.3) are of the size one would expect for an MSE.

**Choice.** MSE $= R^{-1}\sum_r(\hat\tau_r-\tau)^2$ for every estimator. When
comparing with the paper we show DIF MSE separately and add a diagnostic column,
$\mathrm{Var}(Y\mid T=t)+\mathrm{Var}(Y\mid T=t')$, estimated from the simulated
samples.

**Reason.** This is the standard definition and it agrees with the CM/GPS columns.
The inconsistency is reported, not copied.

---

## D16. Interval method for coverage — deferred to Stage 4

**Issue.** Tables 1 and 3 report coverage for CM and GPS (and none for DIF)
without saying how intervals are built. The only hint is Table 5 (Section 6),
where intervals are estimate ± 1.96 × "sd", with "sd" described as "the sample
standard deviation". GPS coverage of 96–100% suggests conservative intervals.

**Choice.** Not specified now. In Stage 4 (optional, after your OK) we will pick
and justify a method derived from the estimator's form, e.g. a variance estimator
derived for 1-NN imputation, or a resampling scheme, and report coverage.

---

## D17. Replications and random numbers — impact L

**Choice.** $R = 500$ per simulation, as in the paper. Random numbers use
`L'Ecuyer-CMRG` with one stream per replication: stream $r$ is
`nextRNGStream` applied $r$ times to a recorded master seed. The draw order within
a replication is fixed: candidate batches (W, then T), then $\varepsilon$. Each
replication's results are saved to its own file as soon as it finishes.

**Reason.** This makes each replication reproducible on its own, safe to run in
parallel on 4 cores, and resumable.

---

## D18. Balance-assessment formulas — out of scope

The balance algorithms (pp. 971, 973) are not used in the simulation tables.
Noted issues, in case they are implemented later as diagnostics:
* Step 1 (p. 971) gives $\bar X_t = \frac1n\sum_{i:T_i=t}X_i$; this should be $1/n_t$.
* Step 2 (p. 971) uses $t$ both as the fixed level and as the summation index in
  $S^2_{X|t} = \frac1A\sum_{t=1}^A S^2_{t,X}$.
* $ND^{GPS}_t$ (p. 973) uses $S_{\Pr(t\mid X,Z)\mid T}$, which is not defined. Step 2
  says "means and sample variance" but gives only the mean. The sign of the
  difference is the opposite of $ND^{cov}_t$.

---

## D19. Measures added beyond the paper

Sim 1 signed bias and MSE (the paper gives |bias| only), empirical SD, and Monte
Carlo SEs of bias and MSE are added for all estimators. The paper reports none
of these.

---

## D20. Unresolved: residual mismatch in the DIF column

Under the chosen DGP (D1–D4) the population DIF bias reproduces the paper's Sim 2
column well (sign-aligned correlation 0.956, slope 1.08, R² 0.91). The RMS
difference, 0.40, is still well above the paper's Monte Carlo error (about
0.02–0.06). In the per-treatment decomposition $h_t$ (spec §3), the differences
are concentrated in $t = 7$ (about 0.6) and $t = 4$ (about 0.4) in Sim 2, and in
$t = 3$ (about 0.7) in Sim 1, using the internally consistent sign pattern. None of
the variants tried (D1, D2) removes them, and Sim 1 and Sim 2 share $\eta_1..\eta_4$
and $\beta_1..\beta_4$ yet disagree on different treatments. So the cause is
unknown: possibly a DGP detail not stated in the paper, an implementation
difference, or covariates held fixed across replications. We do not adjust the
DGP to fit these numbers. Stage 3 will show whether the simulated DIF agrees with
the population calculation, and the remaining gap will be reported as a
discrepancy with the paper.

**Stage 3 status.** At R = 500 the simulated DIF bias agrees with the population
calculation: |z| ≤ 2.4 over the 6 Sim 1 pairs and ≤ 1.9 over the 36 Sim 2 pairs,
in Monte Carlo SE units. Our implementation therefore matches the specified DGP,
and the gap with the paper's DIF column (sign-aligned correlation 0.958, 8 of 36
pairs within 2 SE) is reported as a discrepancy with the paper.

---

## D21. Sim 1 outcomes are not centred — impact L

**Issue.** Sim 1 uses $Y(t) = [x,z]^\top\beta_t$ with raw $Z_2, Z_3$. Sim 2 centres
them.

**Choice.** Take each as written. Sim 1's truth is nonzero (D4).

**Reason.** Adding a constant $c_t$ to $Y(t)$ shifts $\hat\tau$ and $\tau$ by the
same amount for all three estimators, so bias, MSE and coverage are unaffected.

---

# Stage 2 implementation choices

## D22. Fitting the multinomial logit — impact L

**Choice.** A hand-written Newton–Raphson (`R/mlogit.R`). Baseline category 1,
an intercept, and covariates standardised internally (fitted probabilities do
not change). Each step is halved if it would lower the log-likelihood. The fit
counts as converged when the largest Newton step is below 1e-8 on the
standardised scale, with at most 100 iterations. Convergence and iteration
counts are stored for every replication.

**Reason.** It is fully transparent and exact. `nnet::multinom` fits by BFGS
with a loose default tolerance. Tests show the fit equals `glm(binomial)` with
two classes, matches `nnet::multinom` with a tight tolerance on Sim 2-type data
(fitted probabilities within 1e-4, log-likelihood at least as high), satisfies
the score equations to 1e-6, and recovers known coefficients. Pilot: every fit
converged in 5–9 iterations.

## D23. Random draws — impact L

**Choice.** $X = U R$ with $U$ an $N\times 4$ matrix of standard normals and $R$
the upper Cholesky factor of $\Sigma_X$, rather than `MASS::mvrnorm`, whose
eigendecomposition can differ across LAPACK builds. Then
$Z_1$ = `runif(N, -3, 3)`, $Z_2$ = `rchisq(N, 1)`, $Z_3$ = `rbinom(N, 1, .5)`.
$T$ comes from one uniform per candidate, by inverse CDF:
$T = \min\{t : u \le \sum_{s\le t} p_s\}$. Candidates are drawn in batches of
1000 with order (X, Z₁, Z₂, Z₃, u) within a batch. After the groups are filled,
$\varepsilon$ is drawn for the kept units, sorted by treatment and then by order
of drawing. The batch size is part of the random-number protocol (`R/config.R`).

**Note.** Pilot replications 1–20 use the same streams as replications 1–20 of
the full run, so they reproduce exactly.

## D24. Tie tolerance — impact L

**Choice.** A donor counts as tied with the nearest one if its distance is
within $10^{-10}\max(d_{min}, 1)$ of the minimum, and tied donors' outcomes are
averaged (D10).

**Reason.** Floating-point rounding can separate distances that are
mathematically equal; the hand-made test case has exact ties of that kind. With
continuous covariates, genuinely different distances this close do not occur in
practice. A unit's distance to itself is exactly 0, because distances are built
from coordinate differences.

## D25. Sensitivity variants of GPS matching (agreed at the Stage-1 review)

**Choice.** Every replication also computes
* `GPS_cs`: the multinomial logit with main effects plus
  $x_1^2,\dots,x_4^2, z_1^2, z_2^2$. This nests every true $\eta_t$, since each
  equation has its own coefficients. $z_3^2 = z_3$ is omitted as collinear.
* `GPS_true`: matching on the true $p_t(W)$.

The matching step is identical to the primary `GPS` (D13). Pilot note: the
correctly specified model gives some fitted probabilities as small as 1e-28 in
Sim 2 (separation in the $Z_2^2$ tail). These are valid matching scores and are
used as they are.

---

# Stage 4: intervals and coverage

## D26. Interval method — impact H for coverage

**Issue.** Tables 1 and 3 report coverage for CM and GPS but never say how the
intervals were built. The only hint is in Section 6 (Table 5): intervals are
"estimate ± 1.96 × sd", with sd "the sample standard deviation".

**Options.**
(a) A normal-approximation interval with an analytic variance derived from the
estimator's own form (below).
(b) A nonparametric bootstrap, resampling units within each treatment group and
recomputing everything, including the Mahalanobis covariances and the logit fit.
(c) Section 6's recipe, read as ± 1.96 × the SD of the unit-level imputed effects.
(d) "Oracle" intervals, ± 1.96 × the Monte Carlo SD of the estimates across replications.

**Choice.** (a) is the method. (d) is reported only as a diagnostic, because it
uses the true sampling SD, which an analyst never has; it shows how much of any
under-coverage is due to bias alone. For DIF the interval is
$\bar Y_t - \bar Y_{t'} \pm 1.96\,(s_t^2/n_t + s_{t'}^2/n_{t'})^{1/2}$.

**Reason.**
* (a) follows directly from the estimator, is cheap, and accounts for both
  sources of sampling variation under the fixed-group-size design (D3).
* (b) would cost B × about 1 s per Sim 2 replication: roughly 1.5–3 h of compute
  on 4 cores for B = 100–200, more than the time limit allows. It is also on weak
  footing for nearest-neighbour matching: a bootstrap sample contains copies of
  the same unit, and those copies become exact matches for each other, which the
  original sample never has.
* (c) gives a prediction interval for a single unit's effect, not a confidence
  interval for the average. It would be many times too wide (coverage ≈ 100%),
  which does not fit the paper's CM coverage of 92% for some pairs.

**Derivation of (a).** Every estimator here is linear in the observed outcomes:
$\hat\tau(t,t') = n^{-1}\{\sum_{j\in G_t} w_j Y_j - \sum_{j\in G_{t'}} w_j Y_j\}$, where
$w_j$ is the number of units, $j$ itself included, whose imputed $Y(t)$ is $Y_j$
(fractional under ties). Let $V$ be what the estimator matches on: all covariates
for CM, the scalar score for GPS. Write $Y_j = m_{T_j}(V_j) + e_j$ with
$m_t(v) = E[Y \mid V = v, T = t]$. Then

$$\hat\tau - \tau = \underbrace{\Big[\tfrac1n\textstyle\sum_i \tau(V_i) - \tau\Big]}_{\text{units drawn}}
+ \underbrace{B_n}_{\text{matching bias}} + \underbrace{\tfrac1n\textstyle\sum_j c_j e_j}_{\text{outcome noise}},$$

with $c_j = \pm w_j$ and $\tau(V_i)$ the effect at unit $i$'s matching variables.
Given the covariates, the noise term has variance
$\frac{1}{n^2}\sum_j w_j^2\sigma^2(V_j)$. Under D3 the units are independent with
fixed group sizes, so the first term has variance
$\frac{1}{n^2}\sum_s n_s \mathrm{Var}\{\tau(V)\mid T=s\}$. The estimates are:
* $\hat\sigma^2_j = (Y_j - Y_l)^2/2$, where $l$ is $j$'s nearest other unit in the
  same group, measured in the estimator's own metric (ties averaged).
* $\widehat{\mathrm{Var}}\{\tau(V)\mid T=s\}$ is the within-group sample variance of
  the unit-level imputed effects $\hat Y_i(t) - \hat Y_i(t')$, minus their average
  noise variance $\sum_j a_{ij}^2\hat\sigma^2_j$ (both levels), floored at 0.

The interval for the population effect uses
$\hat V = \hat V_{noise} + \hat V_{units}$. The interval for the sample-average
effect $n^{-1}\sum_i\{\mu_t(W_i)-\mu_{t'}(W_i)\}$ uses $\hat V_{noise}$ only.

**Limitations, stated in advance.**
1. The matching bias $B_n$ is not in the interval, by design. The interval
   describes the estimator's sampling variability, so bias lowers coverage. The
   oracle diagnostic (d) separates the two effects.
2. Uncertainty from estimating the propensity score is ignored.
3. For GPS, the noise term and the units term can be slightly correlated, because
   $Y(t)$ and $Y(t')$ are matched on different scores. This is ignored.

**Checks.** The variance pieces reproduce a hand-worked 6-unit example, and in a
randomised setting with a heterogeneous effect the estimated variance falls
within 0.7–1.45 of the Monte Carlo variance for both CM and GPS
(`tests/testthat/test-variance.R`).

**What is reported.** Coverage of the nominal 95% interval against
$\tau_{super}$ (the paper's truth, used for the comparison), against $\tau_{mix}$,
and against the sample-average effect (conditional interval); oracle coverage;
and the SE calibration ratio, mean estimated SE / Monte Carlo SD. With 500
replications, the Monte Carlo SE of a coverage near 95% is about 1 percentage
point.
