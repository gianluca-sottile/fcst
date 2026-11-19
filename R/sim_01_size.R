################################################################################
#  sim_01_size.R
#
#  Empirical size (type I error) of FCST under the null hypothesis of random
#  labelling (beta = 0: type probabilities do not depend on location).
#
#  Factors
#    process : "hpp", "ipp", "thomas"     (ground point process)
#    n       : 200, 500                    (number of points)
#    labels  : "oracle" (true types), "est_easy" (spectrum + FPCA + k-means,
#              SNR = 1, typing almost perfect) or "est_hard" (same typing,
#              SNR = 0.05, ARI about 0.67)
#    weights : "none", "inverse_intensity"
#  Within each (process, n, labels) the two weighting schemes are applied to
#  the same data (paired comparison); all other settings are independent.
#
#  Recorded for each replicate
#    p_global      global test (sum over types, max over radii)
#    p_maxT        first version of the global test (max over types and radii)
#    p_type1_raw   type-specific test for type 1
#    any_holm      1 if at least one type is rejected after Holm (FWER)
#    p_centroid    centroid statistic, min over types with Bonferroni
#
#  Output: results/sim_01_size_raw.csv, results/sim_01_size_summary.csv
################################################################################

source("fcst_methods.R")
dir.create("results", showWarnings = FALSE)

## ---- settings ----------------------------------------------------------------
n_rep     <- 1000           # replicates per setting (paper: 1000)
B         <- 999            # permutations per test  (paper: 999)
n_cores   <- 4              # > 1 runs replicates in parallel (macOS/Linux)
K         <- 3
radii     <- c(0.05, 0.10, 0.15, 0.20)
alpha     <- 0.05
processes <- c("hpp", "ipp", "thomas")
sizes     <- c(200, 500)
label_modes  <- c("oracle", "est_easy", "est_hard")
snr_of_mode  <- c(oracle = NA, est_easy = 1, est_hard = 0.05)
weight_modes <- c("none", "inverse_intensity")
seed_base <- 1e7            # different for every script

## ---- one replicate: both weighting schemes on the same data ----------------------
run_one <- function(process, n, label_mode) {
  X <- sim_points(n, process)
  types <- sim_types(type_probabilities(X, K = K, beta = 0))
  if (label_mode == "oracle") {
    labels <- types
  } else {
    Y <- sim_curves(types, snr = snr_of_mode[[label_mode]])
    labels <- type_curves(Y, K, representation = "spectrum",
                          algorithm = "fpca_kmeans")
  }
  ari <- adjusted_rand_index(types, labels)
  rows <- list()
  for (w in weight_modes) {
    test <- fcst_test(X, labels, radii, B = B, weights = w,
                      centroid = (w == "none"))
    p_centroid <- NA_real_
    if (w == "none") p_centroid <- min(test$centroid$p_bonf, na.rm = TRUE)
    rows[[w]] <- data.frame(process = process, n = n, labels = label_mode,
                            weights = w, ari = ari,
                            p_global = test$p_global,
                            p_maxT = min(test$table$p_maxT, na.rm = TRUE),
                            p_type1_raw = test$table$p_raw[1],
                            any_holm = 1 * any(test$table$p_holm <= alpha,
                                               na.rm = TRUE),
                            p_centroid = p_centroid)
  }
  return(do.call(rbind, rows))
}

## ---- loop over settings -------------------------------------------------------
all_rows <- list()
setting_id <- 0
for (process in processes) {
  for (n in sizes) {
    for (label_mode in label_modes) {
      setting_id <- setting_id + 1
      cat("setting", setting_id, ": process =", process, " n =", n,
          " labels =", label_mode, "\n")
      rows <- run_replicates(function() run_one(process, n, label_mode),
                             n_rep, setting_id, n_cores, seed_base)
      all_rows <- c(all_rows, rows)
    }
  }
}
raw <- do.call(rbind, all_rows)
write.csv(raw, "results/sim_01_size_raw.csv", row.names = FALSE)

## ---- summary: rejection rates at level alpha -------------------------------------
mc_se <- function(p_hat, m) sqrt(p_hat * (1 - p_hat) / m)
settings <- unique(raw[, c("process", "n", "labels", "weights")])
summary_rows <- list()
for (s in 1:nrow(settings)) {
  sel <- raw$process == settings$process[s] & raw$n == settings$n[s] &
    raw$labels == settings$labels[s] & raw$weights == settings$weights[s]
  m <- sum(sel)
  rej_global <- mean(raw$p_global[sel] <= alpha)
  summary_rows[[s]] <- data.frame(
    settings[s, ], n_rep = m,
    mean_ari = mean(raw$ari[sel]),
    size_global = rej_global,
    se_global = mc_se(rej_global, m),
    size_maxT = mean(raw$p_maxT[sel] <= alpha),
    size_type1 = mean(raw$p_type1_raw[sel] <= alpha),
    fwer_holm = mean(raw$any_holm[sel]),
    size_centroid_bonf = mean(raw$p_centroid[sel] <= alpha))
}
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, "results/sim_01_size_summary.csv", row.names = FALSE)
print(summary_table, digits = 3)
## Expected: size_global, size_maxT, size_type1 within alpha +/- 2 se;
## fwer_holm and size_centroid_bonf at or below alpha.
## size_centroid_bonf is NA for inverse weights (not computed).
