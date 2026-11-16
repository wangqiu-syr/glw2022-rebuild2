# The paper's simulation tables, transcribed

Values are exactly as printed in Guo, Liu & Wang (2022). Coverage is in percent,
without the "%" sign.

| File | Source | Rows | Columns |
|---|---|---|---|
| `table1_sim1.csv` | Table 1, p. 974 (Simulation 1) | 6 pairs | `dif_absbias, cm_absbias, cm_coverage, gps_absbias, gps_coverage` |
| `table3_sim2.csv` | Table 3, pp. 976–977 (Simulation 2) | 36 pairs | `dif_bias, dif_mse, cm_bias, cm_mse, cm_coverage, gps_bias, gps_mse, gps_coverage` |

`t`, `tprime` give the pair in $\hat\tau(t, t')$. Table 1 reports absolute bias
only. The signs in Table 3 follow the paper; see `docs/decisions.md` (D5) for why
they appear reversed relative to the paper's own definition.

How the transcription was checked (details in `docs/spec.md` §7):
1. Values were parsed from the PDF text layer and compared cell by cell with
   300-dpi renders of the pages.
2. `scripts/stage1_check_tables.R` checks that every bias column has the form
   $h_t - h_{t'}$, which the estimators imply. The maximum residual is 7e-5,
   i.e. rounding error.
