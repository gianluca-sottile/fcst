################################################################################
#  sim_09_attribution.R
#
#  Attribution of a departure to individual types: FCST against the standard
#  type-specific tools for multitype point patterns (Section 6.6 and Section S8
#  of the Supplementary Material).
#
#  Type-specific tests compared (each followed by Holm's adjustment over the
#  K types, all calibrated with the same kind of random-labelling permutations):
#    fcst      T_k = max_r Z_k(r), radii 0.05, 0.10, 0.15, 0.20 (the paper's test)
#    envelope  global envelope test (Myllymaki et al. 2017) on the type-k pair
#              count U_k(r), equivalent to the type-k K-function K_kk(r) under
#              random labelling, on a grid of 20 radii up to 0.20; extreme rank
#              length (ERL) ordering, one-sided ("greater"), package GET.
#              Without GET, a studentised maximum-deviation envelope is used.
#    nn        nearest-neighbour test for class k (Dixon 2002): number of type-k
#              events whose nearest neighbour is also of type k, calibrated by
#              permutation.
#
#  Design: homogeneous Poisson, n = 300, K = 3; one hotspot, two hotspots, band;
#  beta in {0, 3, 6}; true types and estimated types at SNR = 0.05.
#  Reported: probability that type 1 (the concentrated type) is declared
#  concentrated, that another type is, and that any type is (FWER under H0).
#
#  Output: results/sim_09_attribution_raw.csv, results/sim_09_attribution_summary.csv
################################################################################

source("fcst_methods.R")
dir.create("results", showWarnings = FALSE)
have_get <- requireNamespace("GET", quietly = TRUE)
if (!have_get) message("package GET not installed: studentised envelope used instead")

## ---- settings ----------------------------------------------------------------------
n_rep     <- 500
B         <- 999
n_cores   <- 1
n         <- 300
K         <- 3
radii     <- c(0.05, 0.10, 0.15, 0.20)
r_grid    <- seq(0.01, 0.20, length.out = 20)
alpha     <- 0.05
patterns  <- c("one_bump", "two_bumps", "band")
betas     <- c(0, 3, 6)
label_modes <- c("oracle", "est_hard")
snr_hard  <- 0.05
seed_base <- 9e7

rename_type1 <- function(labels_raw, types) {
  # for REPORTING only: the estimated type overlapping most with true type 1
  # is called type 1 (the tests do not depend on the numbering)
  overlap <- numeric(K)
  for (k in 1:K) overlap[k] <- sum(labels_raw == k & types == 1)
  target <- which.max(overlap)
  labels <- labels_raw
  labels[labels_raw == target] <- 1
  labels[labels_raw == 1] <- target
  return(labels)
}

envelope_p <- function(obs, sims) {
  # obs: vector over r; sims: matrix (number of radii) x B; one-sided "greater"
  if (have_get) {
    cs <- GET::create_curve_set(list(r = r_grid, obs = obs, sim_m = sims))
    res <- GET::global_envelope_test(cs, type = "erl", alternative = "greater")
    return(attr(res, "p"))
  }
  all_curves <- cbind(obs, sims)
  mu <- rowMeans(all_curves)
  s  <- apply(all_curves, 1, sd)
  s[s == 0] <- 1
  stat <- apply((all_curves - mu) / s, 2, max)
  return((1 + sum(stat[-1] >= stat[1])) / (length(stat)))
}

## ---- one replicate ---------------------------------------------------------------
run_one <- function(pattern, beta, label_mode) {
  X <- sim_points(n, "hpp")
  pat <- if (beta == 0) "none" else pattern
  types <- sim_types(type_probabilities(X, K = K, beta = beta, pattern = pat))
  if (label_mode == "oracle") {
    labels <- types
  } else {
    Y <- sim_curves(types, snr = snr_hard)
    labels <- rename_type1(type_curves(Y, K, representation = "spectrum",
                                       algorithm = "fpca_kmeans"), types)
  }

  # FCST
  t0 <- proc.time()["elapsed"]
  p_fcst <- fcst_test(X, labels, radii, B = B)$table$p_holm
  time_fcst <- proc.time()["elapsed"] - t0

  # envelope and nearest-neighbour tests share the same permutations
  t0 <- proc.time()["elapsed"]
  pairs_grid <- build_pair_list(X, r_grid, kernel = "ball")
  D <- as.matrix(dist(X)); diag(D) <- Inf
  nn <- apply(D, 1, which.min)
  nn_count <- function(lab) sapply(1:K, function(k) sum(lab == k & lab[nn] == k))
  U_obs <- pair_statistic(labels, pairs_grid, K)
  N_obs <- nn_count(labels)
  U_null <- array(NA_real_, dim = c(K, length(r_grid), B))
  N_null <- matrix(NA_real_, nrow = B, ncol = K)
  for (b in 1:B) {
    lab_b <- labels[sample.int(n)]
    U_null[, , b] <- pair_statistic(lab_b, pairs_grid, K)
    N_null[b, ] <- nn_count(lab_b)
  }
  p_env <- sapply(1:K, function(k) envelope_p(U_obs[k, ], U_null[k, , ]))
  time_env <- proc.time()["elapsed"] - t0
  p_nn <- sapply(1:K, function(k) (1 + sum(N_null[, k] >= N_obs[k])) / (B + 1))
  p_env <- holm_adjust(p_env)
  p_nn  <- holm_adjust(p_nn)

  out <- data.frame(pattern = pattern, beta = beta, labels = label_mode,
                    ari = adjusted_rand_index(types, labels),
                    time_fcst = as.numeric(time_fcst),
                    time_envelope_nn = as.numeric(time_env))
  for (m in c("fcst", "env", "nn")) {
    p <- switch(m, fcst = p_fcst, env = p_env, nn = p_nn)
    out[[paste0(m, "_type1")]] <- 1 * (p[1] <= alpha)
    out[[paste0(m, "_other")]] <- 1 * any(p[-1] <= alpha)
    out[[paste0(m, "_any")]]   <- 1 * any(p <= alpha)
  }
  return(out)
}

## ---- loop over settings ------------------------------------------------------------
all_rows <- list()
setting_id <- 0
for (pattern in patterns) {
  for (beta in betas) {
    if (beta == 0 && pattern != patterns[1]) next      # the null is run once
    for (label_mode in label_modes) {
      setting_id <- setting_id + 1
      cat("setting", setting_id, ":", ifelse(beta == 0, "null", pattern),
          "beta =", beta, label_mode, "\n")
      rows <- run_replicates(function() run_one(pattern, beta, label_mode),
                             n_rep, setting_id, n_cores, seed_base)
      all_rows <- c(all_rows, rows)
    }
  }
}
raw <- do.call(rbind, all_rows)
raw$pattern[raw$beta == 0] <- "none"
write.csv(raw, "results/sim_09_attribution_raw.csv", row.names = FALSE)

## ---- summary -------------------------------------------------------------------------
keys <- unique(raw[, c("pattern", "beta", "labels")])
rows <- list()
for (s in 1:nrow(keys)) {
  sel <- raw$pattern == keys$pattern[s] & raw$beta == keys$beta[s] &
    raw$labels == keys$labels[s]
  r <- data.frame(keys[s, ], n_rep = sum(sel), mean_ari = mean(raw$ari[sel]))
  for (v in grep("_(type1|other|any)$", names(raw), value = TRUE)) r[[v]] <- mean(raw[[v]][sel])
  r$time_fcst <- median(raw$time_fcst[sel])
  r$time_envelope_nn <- median(raw$time_envelope_nn[sel])
  rows[[s]] <- r
}
summary_table <- do.call(rbind, rows)
write.csv(summary_table, "results/sim_09_attribution_summary.csv", row.names = FALSE)
print(summary_table, digits = 3)
