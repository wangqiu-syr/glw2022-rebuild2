# CM distance-metric sensitivity (D9)

Same 500 + 500 data sets as Stage 3 (same random streams). The primary CM estimates reproduce Stage 3 exactly
(max |difference|: Sim 1 4.9e-15, Sim 2 5.1e-15).

Bias against tau_super; paper's Sim 2 signs multiplied by -1 (D5). Paper CM: Sim 1 mean |bias| 0.1432; Sim 2 mean |bias| 0.1394, mean MSE 0.3677.

| metric | S1 mean abs bias | S1 cor with paper | S1 mean MSE | S2 mean abs bias | S2 cor with paper | S2 slope paper~ours | S2 mean MSE | S2 mean abs bias vs tau_mix |
|---|---|---|---|---|---|---|---|---|
| blocks: Mahalanobis(X) + Mahalanobis(Z) [primary] | 0.1663 | 0.0101 | 0.3644 | 0.2986 | 0.7336 | 0.3247 | 0.2936 | 0.2657 |
| stacked: Mahalanobis(X, Z) | 0.1787 | 0.0742 | 0.3675 | 0.3037 | 0.7573 | 0.3272 | 0.2897 | 0.2702 |
| Mahalanobis(X) + Euclidean(Z) | 0.1460 | -0.1400 | 0.3384 | 0.2657 | 0.7041 | 0.3542 | 0.2919 | 0.2358 |

Per-replication estimates: `estimates_sim1.csv.gz`, `estimates_sim2.csv.gz`. Per-pair values: `performance_sim1.csv`, `performance_sim2.csv`. Figure: `figures/cm_metric_sensitivity_sim2.png`.
