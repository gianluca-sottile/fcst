################################################################################
#  sim_05b_segregation_cv.R
#
#  Fairness check for the comparison of Section 6.6 (Table S2 of the
#  Supplementary Material): does the segregation test gain power when its bandwidth
#  is chosen by cross-validation, as in Diggle et al. (2005) and spatstat,
#  instead of the fixed normal-reference rule used in sim_05?
#
#  The data are generated with the same seeds and the same order of random
#  draws as sim_05. NOTE (11 Nov 2026): in the runs actually performed the
#  replicates do NOT coincide with those of sim_05 at SNR = 0.05 (ARI and
#  p_fcst differ), probably because the two scripts were run in different R
#  sessions/settings. The comparison within this script is unaffected: all
#  tests here are applied to the same data.
#
#  Tests on the estimated types
#    fcst      FCST default
#    seg_fixed segregation test, normal-reference bandwidth (as in sim_05)
#    seg_cv    bandwidth by likelihood cross-validation on the observed labels,
#              then fixed across permutations (spatstat default behaviour)
#    seg_cvp   cross-validation repeated for every permutation (exact, slower)
#
#  Output: results/sim_05b_segregation_cv_raw.csv
#          results/sim_05b_segregation_cv_summary.csv
################################################################################

source("fcst_methods.R")
source("fcst_competitors.R")
dir.create("results", showWarnings = FALSE)

## ---- settings (as in sim_05) ---------------------------------------------------
n_rep      <- 300            # paper: 300
n_rep_null <- 1000           # paper: 1000
B          <- 999
n_cores    <- 1
n          <- 300
K          <- 3
radii      <- c(0.05, 0.10, 0.15, 0.20)
alpha      <- 0.05
patterns   <- c("one_bump", "two_bumps", "band")
betas      <- c(2, 4, 6)
snr_levels <- c(1, 0.05)
seed_base  <- 5e7            # same as sim_05: same data

run_one <- function(pattern, beta, snr) {
  # data generation: identical sequence of random draws to sim_05
  X <- sim_points(n, "hpp")
  types <- sim_types(type_probabilities(X, K = K, beta = beta,
                                        pattern = pattern))
  Y <- sim_curves(types, snr = snr)
  labels <- type_curves(Y, K, representation = "spectrum",
                        algorithm = "fpca_kmeans")
  t_fcst <- fcst_test(X, labels, radii, B = B)
  t_fix  <- segregation_test(X, labels, B = B)
  t_cv   <- segregation_test(X, labels, B = B, bw = "cv")
  t_cvp  <- segregation_test(X, labels, B = B, bw = "cv_perm")
  data.frame(pattern = pattern, beta = beta, snr = snr,
             ari = adjusted_rand_index(types, labels),
             bw_fixed = default_bandwidth(X), bw_cv = t_cv$bw,
             p_fcst = t_fcst$p_global,
             p_seg_fixed = t_fix$p_value,
             p_seg_cv = t_cv$p_value,
             p_seg_cvp = t_cvp$p_value)
}

## ---- loop over settings (same order and numbering as sim_05) ---------------------
all_rows <- list()
setting_id <- 0
for (snr in snr_levels) {
  setting_id <- setting_id + 1
  cat("setting", setting_id, ": null, snr =", snr, "\n")
  rows <- run_replicates(function() run_one("none", 0, snr), n_rep_null,
                         setting_id, n_cores, seed_base)
  all_rows <- c(all_rows, rows)
}
for (pattern in patterns) {
  for (beta in betas) {
    for (snr in snr_levels) {
      setting_id <- setting_id + 1
      cat("setting", setting_id, ":", pattern, "beta =", beta, "snr =", snr, "\n")
      rows <- run_replicates(function() run_one(pattern, beta, snr), n_rep,
                             setting_id, n_cores, seed_base)
      all_rows <- c(all_rows, rows)
    }
  }
}
raw <- do.call(rbind, all_rows)
write.csv(raw, "results/sim_05b_segregation_cv_raw.csv", row.names = FALSE)

## ---- summary ---------------------------------------------------------------------
settings <- unique(raw[, c("pattern", "beta", "snr")])
summary_rows <- list()
for (s in 1:nrow(settings)) {
  sel <- raw$pattern == settings$pattern[s] & raw$beta == settings$beta[s] &
    raw$snr == settings$snr[s]
  summary_rows[[s]] <- data.frame(
    settings[s, ], n_rep = sum(sel), mean_ari = mean(raw$ari[sel]),
    median_bw_fixed = median(raw$bw_fixed[sel]),
    median_bw_cv = median(raw$bw_cv[sel]),
    fcst = mean(raw$p_fcst[sel] <= alpha),
    seg_fixed = mean(raw$p_seg_fixed[sel] <= alpha),
    seg_cv = mean(raw$p_seg_cv[sel] <= alpha),
    seg_cvp = mean(raw$p_seg_cvp[sel] <= alpha))
}
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, "results/sim_05b_segregation_cv_summary.csv",
          row.names = FALSE)
print(summary_table, digits = 3)
