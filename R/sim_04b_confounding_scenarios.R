################################################################################
#  sim_04b_confounding_scenarios.R
#
#  Confounding by a spatially structured covariate: five scenarios
#  (Section 6.7 of the paper).
#
#  The covariate z_i is the distance of event i from a "station" at
#  (0.5, 0.5), the analogue of source-receiver distance in seismology.
#  Aims: (i) a covariate acting on the curves can produce rejections without
#  any shape-location dependence; (ii) stratification on z and invariant
#  representations remove them; (iii) they do not remove genuine
#  shape-location dependence that z does not explain.
#
#  Type model: P(type = k | x) proportional to exp(eta_k(x)), eta_2 = eta_3 = 0,
#  eta_1 given below; g(z) = 1 - z / 0.4 (positive near the station).
#
#  Scenario                          eta_1                       curves
#  1 amplitude      z on amplitude   0                           amp = exp(-4 z)
#  2 type_z         type on z        beta_z g(z)                 standard
#  3 type_z_hot     z + hotspot      beta_z g(z) + beta_h h(x)   standard
#  4 type_z_x       z + abscissa     beta_z g(z) + beta_x (x - 0.5) / 0.29
#  5 snr_z          z on SNR only    0                           SNR = snr0 exp(-c z)
#
#  Expected (spectrum typing):
#                 unrestricted     stratified on z
#  1 amplitude    level            level      (raw curves: unrestricted rejects)
#  2 type_z       rejects          level      (up to strata resolution)
#  3 type_z_hot   rejects          rejects    (hotspot not explained by z)
#  4 type_z_x     rejects          rejects    (correct: x not used in strata)
#  5 snr_z        may reject       level      (differential typing errors)
#
#  Analyses: spectrum typing with unrestricted permutations and with
#  permutations within 8 and 16 quantile strata of z; in scenario 1 also
#  the raw-curve typing (as in sim_04).
#
#  Output: results/sim_04b_confounding_raw.csv
#          results/sim_04b_confounding_summary.csv
################################################################################

source("fcst_methods.R")
dir.create("results", showWarnings = FALSE)

## ---- settings ----------------------------------------------------------------
n_rep_null <- 1000           # scenarios 1, 2, 5 (no residual dependence)
n_rep_alt  <- 500            # scenarios 3, 4
B        <- 999
n_cores  <- 1
n        <- 300
K        <- 3
radii    <- c(0.05, 0.10, 0.15, 0.20)
alpha    <- 0.05
gamma_z  <- 4                # amplitude decay (scenario 1, as in sim_04)
beta_z   <- 2                # effect of z on type 1 (scenarios 2-4)
beta_h   <- 4                # hotspot not explained by z (scenario 3)
beta_x   <- 1                # second covariate, not stratified (scenario 4)
snr0     <- 1                # SNR at the station (scenario 5)
snr_decay <- 6               # SNR = snr0 * exp(-snr_decay * z): about 0.015 far away
scenarios <- c("amplitude", "type_z", "type_z_hot", "type_z_x", "snr_z")
seed_base <- 4.5e7

g_of_z <- function(z) 1 - z / 0.4

type_prob_scenario <- function(X, z, scenario) {
  eta <- matrix(0, nrow = nrow(X), ncol = K)
  if (scenario == "type_z") {
    eta[, 1] <- beta_z * g_of_z(z)
  } else if (scenario == "type_z_hot") {
    eta[, 1] <- beta_z * g_of_z(z) + beta_h * hotspot_surface(X, "one_bump")
  } else if (scenario == "type_z_x") {
    eta[, 1] <- beta_z * g_of_z(z) + beta_x * (X[, 1] - 0.5) / 0.29
  }
  P <- exp(eta)
  return(P / rowSums(P))
}

## ---- one replicate ---------------------------------------------------------------
run_one <- function(scenario) {
  X <- sim_points(n, "hpp")
  z <- sqrt((X[, 1] - 0.5)^2 + (X[, 2] - 0.5)^2)
  types <- sim_types(type_prob_scenario(X, z, scenario))
  amp <- rep(1, n)
  snr <- 1
  if (scenario == "amplitude") amp <- exp(-gamma_z * z)
  if (scenario == "snr_z") snr <- snr0 * exp(-snr_decay * z)
  Y <- sim_curves(types, snr = snr, amp_multiplier = amp)

  lab_spec <- type_curves(Y, K, representation = "spectrum",
                          algorithm = "fpca_kmeans")
  strata8  <- make_strata(z, n_bins = 8)
  strata16 <- make_strata(z, n_bins = 16)
  near <- z <= median(z)
  row <- data.frame(
    scenario = scenario,
    ari_spec = adjusted_rand_index(types, lab_spec),
    ari_near = adjusted_rand_index(types[near], lab_spec[near]),
    ari_far  = adjusted_rand_index(types[!near], lab_spec[!near]),
    p_spec_unrestricted = fcst_test(X, lab_spec, radii, B = B)$p_global,
    p_spec_strata8  = fcst_test(X, lab_spec, radii, B = B, strata = strata8)$p_global,
    p_spec_strata16 = fcst_test(X, lab_spec, radii, B = B, strata = strata16)$p_global,
    ari_raw = NA, p_raw_unrestricted = NA, p_raw_strata8 = NA,
    p_raw_strata16 = NA)
  if (scenario == "amplitude") {
    lab_raw <- type_curves(Y, K, representation = "raw",
                           algorithm = "fpca_kmeans")
    row$ari_raw <- adjusted_rand_index(types, lab_raw)
    row$p_raw_unrestricted <- fcst_test(X, lab_raw, radii, B = B)$p_global
    row$p_raw_strata8  <- fcst_test(X, lab_raw, radii, B = B, strata = strata8)$p_global
    row$p_raw_strata16 <- fcst_test(X, lab_raw, radii, B = B, strata = strata16)$p_global
  }
  return(row)
}

## ---- loop over scenarios -----------------------------------------------------------
all_rows <- list()
for (s in seq_along(scenarios)) {
  scenario <- scenarios[s]
  n_rep <- n_rep_null
  if (scenario %in% c("type_z_hot", "type_z_x")) n_rep <- n_rep_alt
  cat("scenario", s, ":", scenario, "(", n_rep, "replicates )\n")
  rows <- run_replicates(function() run_one(scenario), n_rep, s, n_cores,
                         seed_base)
  all_rows <- c(all_rows, rows)
}
raw <- do.call(rbind, all_rows)
write.csv(raw, "results/sim_04b_confounding_raw.csv", row.names = FALSE)

## ---- summary ---------------------------------------------------------------------
reject <- function(p) mean(p <= alpha)
summary_rows <- list()
for (scenario in scenarios) {
  sel <- raw$scenario == scenario
  summary_rows[[scenario]] <- data.frame(
    scenario = scenario, n_rep = sum(sel),
    ari_spec = mean(raw$ari_spec[sel]),
    ari_near = mean(raw$ari_near[sel]),
    ari_far  = mean(raw$ari_far[sel]),
    spec_unrestricted = reject(raw$p_spec_unrestricted[sel]),
    spec_strata8      = reject(raw$p_spec_strata8[sel]),
    spec_strata16     = reject(raw$p_spec_strata16[sel]),
    ari_raw           = mean(raw$ari_raw[sel]),
    raw_unrestricted  = reject(raw$p_raw_unrestricted[sel]),
    raw_strata8       = reject(raw$p_raw_strata8[sel]),
    raw_strata16      = reject(raw$p_raw_strata16[sel]))
}
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, "results/sim_04b_confounding_summary.csv", row.names = FALSE)
print(summary_table, digits = 3)
