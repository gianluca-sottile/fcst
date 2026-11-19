################################################################################
#  sim_08_robustness.R
#
#  Robustness of power to the ground process and to the hotspot profile
#  (Section S6 of the Supplementary Material).
#
#  Ground processes (n = 300 points)
#    hpp             homogeneous Poisson
#    ipp             inhomogeneous Poisson, intensity proportional to 0.1 + x
#    thomas          Thomas process, 15 parents placed uniformly
#    thomas_overlap  Thomas process, one parent at each hotspot centre of
#                    type 1 (+ uniform parents up to 15): clusters of events
#                    and concentration of the type overlap
#  Hotspot profile of type 1
#    gauss  Gaussian bump (as in sim_02)
#    t      bivariate Student t (df = 3): same peak and same integral,
#           denser core and heavier tails ("denser hotspots")
#  Geometry: one bump, two bumps.  beta in {0, 3, 6}.
#  Typings: true types (oracle) and estimated types at SNR = 0.05
#  (Fourier modulus + FPCA + k-means, ARI about 0.67).
#  Pair weights: unit and inverse-intensity, applied to the SAME data.
#
#  Under beta = 0 the profile plays no role, so the null is run once per
#  (process, geometry, typing); the geometry still matters for
#  "thomas_overlap", whose parents sit at the hotspot centres.
#
#  Output: results/sim_08_robustness_raw.csv
#          results/sim_08_robustness_summary.csv
################################################################################

source("fcst_methods.R")
dir.create("results", showWarnings = FALSE)

## ---- settings ----------------------------------------------------------------
n_rep     <- 500             # paper: 500
B         <- 999             # paper: 999
n_cores   <- 1
n         <- 300
K         <- 3
radii     <- c(0.05, 0.10, 0.15, 0.20)
alpha     <- 0.05
processes <- c("hpp", "ipp", "thomas", "thomas_overlap")
geometries <- c("one_bump", "two_bumps")
profiles  <- c("gauss", "t")
betas     <- c(3, 6)         # beta = 0 handled separately
label_modes <- c("oracle", "est_hard")
snr_hard  <- 0.05
seed_base <- 8e7

## ---- helpers -------------------------------------------------------------------
pattern_name <- function(geometry, profile) {
  # "one_bump" + "t" -> "one_bump_t"; Gaussian profile keeps the old name
  if (profile == "t") return(paste0(geometry, "_t"))
  return(geometry)
}

simulate_ground <- function(process, geometry) {
  if (process == "thomas_overlap") {
    return(sim_points_thomas(n, n_parents = 15, sigma = 0.03,
                             fixed_parents = hotspot_centres(geometry)))
  }
  return(sim_points(n, process))
}

rename_type1 <- function(labels_raw, types) {
  # For REPORTING only: the estimated type overlapping most with true type 1
  # is renamed 1 (global p-values do not depend on the renaming).
  overlap <- numeric(K)
  for (k in 1:K) overlap[k] <- sum(labels_raw == k & types == 1)
  target <- which.max(overlap)
  labels <- labels_raw
  labels[labels_raw == target] <- 1
  labels[labels_raw == 1] <- target
  return(labels)
}

## ---- one replicate ---------------------------------------------------------------
run_one <- function(process, geometry, profile, beta, label_mode) {
  X <- simulate_ground(process, geometry)
  pattern <- pattern_name(geometry, profile)
  if (beta == 0) pattern <- "none"
  types <- sim_types(type_probabilities(X, K = K, beta = beta,
                                        pattern = pattern))
  if (label_mode == "oracle") {
    labels <- types
  } else {
    Y <- sim_curves(types, snr = snr_hard)
    labels_raw <- type_curves(Y, K, representation = "spectrum",
                              algorithm = "fpca_kmeans")
    labels <- rename_type1(labels_raw, types)
  }
  rows <- list()
  for (w in c("none", "inverse_intensity")) {
    test <- fcst_test(X, labels, radii, B = B, weights = w)
    holm <- test$table$p_holm
    rows[[w]] <- data.frame(process = process, geometry = geometry,
                            profile = profile, beta = beta,
                            labels = label_mode, weights = w,
                            ari = adjusted_rand_index(types, labels),
                            p_global = test$p_global,
                            holm_type1 = 1 * (!is.na(holm[1]) && holm[1] <= alpha),
                            holm_other = 1 * any(holm[-1] <= alpha, na.rm = TRUE))
  }
  return(do.call(rbind, rows))
}

## ---- loop over settings ------------------------------------------------------------
all_rows <- list()
setting_id <- 0
for (process in processes) {
  for (geometry in geometries) {
    for (label_mode in label_modes) {
      # null hypothesis (profile irrelevant; recorded as "none")
      setting_id <- setting_id + 1
      cat("setting", setting_id, ":", process, geometry, "beta = 0", label_mode, "\n")
      rows <- run_replicates(function() run_one(process, geometry, "none", 0,
                                                label_mode),
                             n_rep, setting_id, n_cores, seed_base)
      all_rows <- c(all_rows, rows)
      for (profile in profiles) {
        for (beta in betas) {
          setting_id <- setting_id + 1
          cat("setting", setting_id, ":", process, geometry, profile,
              "beta =", beta, label_mode, "\n")
          rows <- run_replicates(function() run_one(process, geometry, profile,
                                                    beta, label_mode),
                                 n_rep, setting_id, n_cores, seed_base)
          all_rows <- c(all_rows, rows)
        }
      }
    }
  }
}
raw <- do.call(rbind, all_rows)
write.csv(raw, "results/sim_08_robustness_raw.csv", row.names = FALSE)

## ---- summary ---------------------------------------------------------------------
keys <- c("process", "geometry", "profile", "beta", "labels", "weights")
settings <- unique(raw[, keys])
summary_rows <- list()
for (s in 1:nrow(settings)) {
  sel <- rep(TRUE, nrow(raw))
  for (key in keys) sel <- sel & raw[[key]] == settings[[key]][s]
  summary_rows[[s]] <- data.frame(
    settings[s, ], n_rep = sum(sel),
    mean_ari     = mean(raw$ari[sel]),
    power_global = mean(raw$p_global[sel] <= alpha),
    holm_type1   = mean(raw$holm_type1[sel]),
    holm_other   = mean(raw$holm_other[sel]))
}
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, "results/sim_08_robustness_summary.csv", row.names = FALSE)
print(summary_table, digits = 3)
