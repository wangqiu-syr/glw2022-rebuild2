# Pilot run (Stage 2)

R 4.3.3; 4 cores; generated 2026-11-16.

## sim1 (20 replications)

Averages over all pairs (bias against tau_super; with only 20 replications these are rough):

| estimator | mean abs bias | mean MSE | mean MC SE of bias |
|---|---|---|---|
| DIF | 0.611 | 1.761 | 0.248 |
| CM | 0.190 | 0.415 | 0.134 |
| GPS | 0.107 | 0.635 | 0.177 |
| GPS_cs | 0.201 | 0.808 | 0.193 |
| GPS_true | 0.065 | 1.540 | 0.273 |

DIF sanity check: simulated DIF bias minus population DIF bias, in MC-SE units: mean -0.24, max |z| 0.45 over 6 pairs.

## sim2 (20 replications)

Averages over all pairs (bias against tau_super; with only 20 replications these are rough):

| estimator | mean abs bias | mean MSE | mean MC SE of bias |
|---|---|---|---|
| DIF | 1.042 | 2.366 | 0.204 |
| CM | 0.317 | 0.282 | 0.084 |
| GPS | 0.148 | 0.599 | 0.168 |
| GPS_cs | 0.285 | 0.970 | 0.200 |
| GPS_true | 0.327 | 1.842 | 0.287 |

DIF sanity check: simulated DIF bias minus population DIF bias, in MC-SE units: mean -0.33, max |z| 2.94 over 36 pairs.

## Timing and convergence

| sim | reps | s/rep (mean) | s/rep (max) | projected 500 reps on 4 cores (min) | candidates drawn (mean) | main logit converged | cs logit converged | max iter main / cs | min p-hat (cs) |
|---|---|---|---|---|---|---|---|---|---|
| sim1 | 20 | 0.13 | 0.18 | 0.3 | 1000 | 20/20 | 20/20 | 6 / 8 | 3.5e-08 |
| sim2 | 20 | 1.15 | 1.26 | 2.4 | 4050 | 20/20 | 20/20 | 5 / 9 | 2.4e-28 |
