## Simulation settings (docs/spec.md sections 2 and 6).

SIMS <- list(
  sim1 = list(name = "sim1", A = 4L, n_per_group = 100L, master_seed = 20220401L),
  sim2 = list(name = "sim2", A = 9L, n_per_group = 200L, master_seed = 20220402L)
)

N_REPS <- 500L

## Estimators computed in every replication.
##   DIF      sample difference                         (spec 3.1)
##   CM       covariates matching                       (spec 3.2)
##   GPS      GPS matching, main-effects multinomial logit (spec 3.3, D12 primary)
##   GPS_cs   GPS matching, correctly specified logit   (D12 sensitivity, agreed for Stage 3)
##   GPS_true GPS matching on the true probabilities    (D12 sensitivity, agreed for Stage 3)
ESTIMATORS <- c("DIF", "CM", "GPS", "GPS_cs", "GPS_true")

## Candidate batch size used when filling the treatment groups (spec 2.4).
## Fixed, because it determines the order of random draws (D17).
BATCH_SIZE <- 1000L
