################################################################################
#  sim_03_misclassification.R
#
#  Effect of typing errors on power (Lemma 1 and Proposition 3 of the paper).
#  True types are generated with a type-1 hotspot; the labels passed to the
#  test are the true labels corrupted by symmetric misclassification with
#  error rate eps. Theory: the spatial signal functional Delta_1 is
#  multiplied exactly by gamma^2, gamma = 1 - eps * K / (K - 1).
#
#  Part A (deterministic check): Delta computed from p_hat equals
#          gamma^2 * Delta computed from p, for every radius.
#  Part B (Monte Carlo): power of the global and type-1 (Holm) tests versus eps.
#
#  Output: results/sim_03_misclass_summary.csv
################################################################################

source("fcst_methods.R")
dir.create("results", showWarnings = FALSE)

## ---- settings ----------------------------------------------------------------
n_rep   <- 500            # paper: 500
B       <- 999            # paper: 999
n_cores <- 4
n       <- 300
K       <- 3
beta    <- 2.5
radii   <- c(0.05, 0.10, 0.15, 0.20)
alpha   <- 0.05
eps_grid <- c(0, 0.1, 0.2, 0.3, 0.4)
seed_base <- 3e7

## ---- Part A: exact attenuation identity ------------------------------------------
set.seed(seed_base)
X <- sim_points(n, "hpp")
P <- type_probabilities(X, K = K, beta = beta, pattern = "one_bump")
pairs <- build_pair_list(X, radii)
cat("Part A: Delta(p_hat) / Delta(p) versus gamma^2, for each radius\n")
for (eps in eps_grid) {
  g <- attenuation_gamma(eps, K)
  p_hat <- eps / (K - 1) + g * P[, 1]
  ratios <- numeric(length(radii))
  for (m in seq_along(radii)) {
    ratios[m] <- delta_signal(p_hat, pairs, m) / delta_signal(P[, 1], pairs, m)
  }
  cat(sprintf("  eps = %.2f   gamma^2 = %.6f   ratios = %s\n", eps, g^2,
              paste(sprintf("%.6f", ratios), collapse = " ")))
}

## ---- Part B: power versus eps --------------------------------------------------
run_one <- function(eps) {
  X <- sim_points(n, "hpp")
  P <- type_probabilities(X, K = K, beta = beta, pattern = "one_bump")
  types <- sim_types(P)
  labels <- add_label_noise(types, K, eps)
  test <- fcst_test(X, labels, radii, B = B)
  data.frame(eps = eps, p_global = test$p_global,
             p_type1_holm = test$table$p_holm[1])
}

all_rows <- list()
for (e in seq_along(eps_grid)) {
  cat("eps =", eps_grid[e], "\n")
  rows <- run_replicates(function() run_one(eps_grid[e]), n_rep, e, n_cores,
                         seed_base)
  all_rows <- c(all_rows, rows)
}
raw <- do.call(rbind, all_rows)

summary_rows <- list()
for (e in seq_along(eps_grid)) {
  sel <- raw$eps == eps_grid[e]
  summary_rows[[e]] <- data.frame(eps = eps_grid[e],
                                  gamma2 = attenuation_gamma(eps_grid[e], K)^2,
                                  power_global = mean(raw$p_global[sel] <= alpha),
                                  power_type1 = mean(raw$p_type1_holm[sel] <= alpha,
                                                     na.rm = TRUE))
}
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, "results/sim_03_misclass_summary.csv", row.names = FALSE)
print(summary_table, digits = 3)
