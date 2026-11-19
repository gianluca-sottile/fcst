################################################################################
#  sim_02_power.R
#
#  Power of FCST against type-specific spatial concentration.
#  Type 1 has probability p_1(x) proportional to exp(beta * h(x)), where h is
#  one bump, two bumps or an elongated band (see hotspot_surface()).
#
#  Compared statistics (same data, same typing)
#    global   FCST default: sum over types, max over radii (ball kernel)
#    gauss    same with the Gaussian kernel (sensitivity analysis)
#    maxT     max over types and radii (first version of the statistic)
#    centroid centroid statistic S(k), min over types with Bonferroni
#  Attribution (type-specific tests, Holm across types)
#    holm_type1  type 1 (the concentrated type) rejected
#    holm_other  type 2 or 3 rejected (their probabilities also vary in
#                space, since they are depleted where type 1 is in excess)
#
#  Output: results/sim_02_power_raw.csv, results/sim_02_power_summary.csv
################################################################################

source("fcst_methods.R")
dir.create("results", showWarnings = FALSE)

## ---- settings ----------------------------------------------------------------
n_rep    <- 500              # paper: 500
B        <- 999              # paper: 999
n_cores  <- 6
n        <- 300
K        <- 3
radii    <- c(0.05, 0.10, 0.15, 0.20)
bandwidths <- c(1, 2, 4) / 32
alpha    <- 0.05
patterns <- c("one_bump", "two_bumps", "band")
betas    <- c(0, 1, 2, 3, 4, 6)
label_modes <- c("oracle", "est_easy", "est_hard", "est_very_hard")
# SNR per mode: ARI about 1, 1, 0.67 and 0.36
snr_of_mode <- c(oracle = NA, est_easy = 1, est_hard = 0.05,
                 est_very_hard = 0.03)
seed_base <- 2e7

run_one <- function(pattern, beta, label_mode) {
  X <- sim_points(n, "hpp")
  P <- type_probabilities(X, K = K, beta = beta, pattern = pattern)
  types <- sim_types(P)
  if (label_mode == "oracle") {
    labels <- types
  } else {
    Y <- sim_curves(types, snr = snr_of_mode[[label_mode]])
    labels_raw <- type_curves(Y, K, representation = "spectrum",
                              algorithm = "fpca_kmeans")
    # For REPORTING the type-specific results only, the estimated type that
    # overlaps most with true type 1 is renamed 1 (the analyst does not know
    # it; global p-values are invariant to this renaming).
    overlap <- numeric(K)
    for (k in 1:K) overlap[k] <- sum(labels_raw == k & types == 1)
    target <- which.max(overlap)
    labels <- labels_raw
    labels[labels_raw == target] <- 1
    labels[labels_raw == 1] <- target
  }
  test  <- fcst_test(X, labels, radii, B = B, kernel = "ball", centroid = TRUE)
  testg <- fcst_test(X, labels, bandwidths, B = B, kernel = "gaussian")
  holm <- test$table$p_holm
  data.frame(pattern = pattern, beta = beta, labels = label_mode,
             ari = adjusted_rand_index(types, labels),
             p_global = test$p_global,
             p_gauss = testg$p_global,
             p_maxT = min(test$table$p_maxT, na.rm = TRUE),
             p_centroid = min(test$centroid$p_bonf, na.rm = TRUE),
             holm_type1 = 1 * (!is.na(holm[1]) && holm[1] <= alpha),
             holm_other = 1 * any(holm[-1] <= alpha, na.rm = TRUE))
}

all_rows <- list()
setting_id <- 0
for (pattern in patterns) {
  for (beta in betas) {
    for (label_mode in label_modes) {
      setting_id <- setting_id + 1
      cat("setting", setting_id, ": pattern =", pattern, " beta =", beta,
          " labels =", label_mode, "\n")
      rows <- run_replicates(function() run_one(pattern, beta, label_mode),
                             n_rep, setting_id, n_cores, seed_base)
      all_rows <- c(all_rows, rows)
    }
  }
}
raw <- do.call(rbind, all_rows)
write.csv(raw, "results/sim_02_power_raw.csv", row.names = FALSE)

## ---- summary -------------------------------------------------------------------
settings <- unique(raw[, c("pattern", "beta", "labels")])
summary_rows <- list()
for (s in 1:nrow(settings)) {
  sel <- raw$pattern == settings$pattern[s] & raw$beta == settings$beta[s] &
    raw$labels == settings$labels[s]
  summary_rows[[s]] <- data.frame(
    settings[s, ], n_rep = sum(sel),
    mean_ari       = mean(raw$ari[sel]),
    power_global   = mean(raw$p_global[sel] <= alpha),
    power_gauss    = mean(raw$p_gauss[sel] <= alpha),
    power_maxT     = mean(raw$p_maxT[sel] <= alpha),
    power_centroid = mean(raw$p_centroid[sel] <= alpha),
    holm_type1     = mean(raw$holm_type1[sel]),
    holm_other     = mean(raw$holm_other[sel]))
}
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, "results/sim_02_power_summary.csv", row.names = FALSE)
print(summary_table, digits = 3)
