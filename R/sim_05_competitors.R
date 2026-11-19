################################################################################
#  sim_05_competitors.R
#
#  Size and power of FCST compared with competing tests, on the same
#  simulated data (Section 6.6 of the paper).
#
#  Data: homogeneous Poisson pattern, n = 300, K = 3 types; type 1
#  concentrated with strength beta in one hotspot, two hotspots or a band.
#  beta = 0 (null hypothesis) does not depend on the pattern: it is run once
#  per SNR with n_rep_null replicates. Curves as in sim_02, at SNR = 1
#  (typing almost perfect) and SNR = 0.05 (ARI about 0.67).
#
#  Tests (all calibrated by random labelling with B permutations)
#    on the estimated type labels:
#      fcst        FCST default (sum over types, max over radii)
#      fcst_maxT   first version (max over types and radii)
#      centroid    centroid statistic, min over types with Bonferroni
#      segregation Diggle, Zheng & Durr (2005)
#      mscan       multinomial scan statistic (Jung et al. 2010)
#    on the curves (Fourier-modulus representation, no typing):
#      fscan       functional scan statistic (FPCA scores)
#      markK       global mark-weighted K-function with functional marks
#      dcov        distance covariance between locations and curves
#  Computing time (seconds, elapsed) of each test is recorded.
#
#  Output: results/sim_05_competitors_raw.csv
#          results/sim_05_competitors_summary.csv
################################################################################

source("fcst_methods.R")
source("fcst_competitors.R")
dir.create("results", showWarnings = FALSE)

## ---- settings ----------------------------------------------------------------
n_rep      <- 300            # paper: 300
n_rep_null <- 1000           # paper: 1000
B          <- 999            # paper: 999
n_cores    <- 4
n          <- 300
K          <- 3
radii      <- c(0.05, 0.10, 0.15, 0.20)
alpha      <- 0.05
max_frac   <- 0.5            # largest scan window: half of the points
patterns   <- c("one_bump", "two_bumps", "band")
betas      <- c(2, 4, 6)     # beta = 0 handled separately
snr_levels <- c(1, 0.05)
seed_base  <- 5e7

elapsed <- function(expr_time) as.numeric(expr_time["elapsed"])

## ---- one replicate ----------------------------------------------------------
run_one <- function(pattern, beta, snr) {
  X <- sim_points(n, "hpp")
  types <- sim_types(type_probabilities(X, K = K, beta = beta,
                                        pattern = pattern))
  Y <- sim_curves(types, snr = snr)
  Z <- represent_curves(Y, method = "spectrum")
  labels <- type_curves(Y, K, representation = "spectrum",
                        algorithm = "fpca_kmeans")
  O <- neighbour_order(X, max_frac)             # shared by the two scans

  tm <- numeric(0)
  tm["fcst"] <- elapsed(system.time(
    t_fcst <- fcst_test(X, labels, radii, B = B)))
  tm["centroid"] <- elapsed(system.time(
    t_cen <- fcst_test(X, labels, radii, B = B, centroid = TRUE)))
  tm["segregation"] <- elapsed(system.time(
    t_seg <- segregation_test(X, labels, B = B)))
  tm["mscan"] <- elapsed(system.time(
    t_msc <- multinomial_scan(X, labels, B = B, O = O)))
  tm["fscan"] <- elapsed(system.time(
    t_fsc <- functional_scan(X, Z, B = B, O = O)))
  tm["markK"] <- elapsed(system.time(
    t_mk <- mark_weighted_K(X, Z, radii, B = B)))
  tm["dcov"] <- elapsed(system.time(
    t_dcov <- dcov_test(X, Z, B = B)))

  row <- data.frame(pattern = pattern, beta = beta, snr = snr,
                    ari = adjusted_rand_index(types, labels),
                    p_fcst = t_fcst$p_global,
                    p_fcst_maxT = min(t_fcst$table$p_maxT, na.rm = TRUE),
                    p_centroid = min(t_cen$centroid$p_bonf, na.rm = TRUE),
                    p_segregation = t_seg$p_value,
                    p_mscan = t_msc$p_value,
                    p_fscan = t_fsc$p_value,
                    p_markK = t_mk$p_value,
                    p_dcov = t_dcov$p_value)
  for (t in names(tm)) row[[paste0("time_", t)]] <- tm[[t]]
  # (time_centroid includes the pair statistic computed in the same call)
  return(row)
}

## ---- loop over settings -------------------------------------------------------
all_rows <- list()
setting_id <- 0
for (snr in snr_levels) {                       # null cells
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
      cat("setting", setting_id, ": pattern =", pattern, " beta =", beta,
          " snr =", snr, "\n")
      rows <- run_replicates(function() run_one(pattern, beta, snr), n_rep,
                             setting_id, n_cores, seed_base)
      all_rows <- c(all_rows, rows)
    }
  }
}
raw <- do.call(rbind, all_rows)
write.csv(raw, "results/sim_05_competitors_raw.csv", row.names = FALSE)

## ---- summary: rejection rates and median computing times -------------------------
tests <- c("fcst", "fcst_maxT", "centroid", "segregation", "mscan", "fscan",
           "markK", "dcov")
timed <- c("fcst", "centroid", "segregation", "mscan", "fscan", "markK", "dcov")
settings <- unique(raw[, c("pattern", "beta", "snr")])
summary_rows <- list()
for (s in 1:nrow(settings)) {
  sel <- raw$pattern == settings$pattern[s] & raw$beta == settings$beta[s] &
    raw$snr == settings$snr[s]
  row <- data.frame(settings[s, ], n_rep = sum(sel),
                    mean_ari = mean(raw$ari[sel]))
  for (t in tests) {
    row[[t]] <- mean(raw[sel, paste0("p_", t)] <= alpha)
  }
  for (t in timed) {
    row[[paste0("time_", t)]] <- median(raw[sel, paste0("time_", t)])
  }
  summary_rows[[s]] <- row
}
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, "results/sim_05_competitors_summary.csv",
          row.names = FALSE)
print(summary_table, digits = 3)
