################################################################################
#  sim_10_semisynthetic.R
#
#  Semi-synthetic experiment with REAL waveforms (Section S10 of the
#  Supplementary Material). Requires results/app_stead_data.rds and
#  results/app_stead_labels.rds (app_01 and app_02).
#
#  Real curves, simulated locations:
#    * n = 300 locations from a homogeneous Poisson process on [0, 1]^2;
#    * each location receives a waveform family k with probabilities
#      q_k(x) proportional to exp{beta h(x)} for family 1 and 1 otherwise
#      (one hotspot, two hotspots or band, as in Section 6);
#    * the curve placed at the location is drawn at random from the STEAD
#      events of that family (families = types of the primary clustering of the
#      application: family 1 = low-frequency type, family 2 = ordinary
#      seismicity, family 3 = compact group of similar events).
#  Under beta = 0 the waveforms are randomly labelled, so the size of the test
#  is checked on real curves. The curves are then clustered again (the analyst
#  does not know the families), with
#    * k-means with K = 3 on the FPCA scores of the Fourier modulus, and
#    * the joint test over k-means with K = 2, ..., 5 (K unknown).
#
#  Output: results/sim_10_semisynthetic_raw.csv, results/sim_10_semisynthetic_summary.csv
################################################################################

source("fcst_methods.R")
dir.create("results", showWarnings = FALSE)

n_rep     <- 500
B         <- 999
n_cores   <- 1
n         <- 300
radii     <- c(0.05, 0.10, 0.15, 0.20)
alpha     <- 0.05
patterns  <- c("one_bump", "two_bumps", "band")
betas     <- c(2, 4, 6)
K_family  <- 2:5
seed_base <- 1e8

stead <- readRDS("results/app_stead_data.rds")
lab   <- readRDS("results/app_stead_labels.rds")
primary <- lab$labels[, lab$primary]
# family 1 = primary type 2 (low frequency, the one that is concentrated under
# the alternative), family 2 = primary type 1, family 3 = primary type 3
family_members <- list(which(primary == 2), which(primary == 1), which(primary == 3))
Y_all <- stead$Y

run_one <- function(pattern, beta) {
  X <- sim_points(n, "hpp")
  pat <- if (beta == 0) "none" else pattern
  fam <- sim_types(type_probabilities(X, K = 3, beta = beta, pattern = pat))
  pick <- integer(n)
  for (k in 1:3) {
    idx <- which(fam == k)
    if (length(idx) > 0) {
      pick[idx] <- family_members[[k]][sample.int(length(family_members[[k]]),
                                                  length(idx), replace = TRUE)]
    }
  }
  Y <- Y_all[pick, , drop = FALSE]
  lab3 <- type_curves(Y, 3, representation = "spectrum", algorithm = "fpca_kmeans")
  t3 <- fcst_test(X, lab3, radii, B = B)
  # type matching family 1 (for reporting only)
  overlap <- sapply(1:3, function(k) sum(lab3 == k & fam == 1))
  k1 <- which.max(overlap)
  labs <- sapply(K_family, function(K) type_curves(Y, K, representation = "spectrum",
                                                    algorithm = "fpca_kmeans"))
  tj <- fcst_test(X, labs, radii, B = B)
  data.frame(pattern = pattern, beta = beta,
             ari = adjusted_rand_index(fam, lab3),
             p_K3 = t3$p_global,
             holm_family1 = 1 * (t3$table$p_holm[k1] <= alpha),
             holm_other = 1 * any(t3$table$p_holm[-k1] <= alpha),
             p_joint_K2to5 = tj$p_global)
}

all_rows <- list()
setting_id <- 1
cat("setting 1: null\n")
all_rows <- c(all_rows, run_replicates(function() run_one("none", 0), 1000, setting_id,
                                       n_cores, seed_base))
for (pattern in patterns) {
  for (beta in betas) {
    setting_id <- setting_id + 1
    cat("setting", setting_id, ":", pattern, "beta =", beta, "\n")
    all_rows <- c(all_rows, run_replicates(function() run_one(pattern, beta), n_rep,
                                           setting_id, n_cores, seed_base))
  }
}
raw <- do.call(rbind, all_rows)
write.csv(raw, "results/sim_10_semisynthetic_raw.csv", row.names = FALSE)

keys <- unique(raw[, c("pattern", "beta")])
rows <- list()
for (s in 1:nrow(keys)) {
  sel <- raw$pattern == keys$pattern[s] & raw$beta == keys$beta[s]
  rows[[s]] <- data.frame(keys[s, ], n_rep = sum(sel), mean_ari = mean(raw$ari[sel]),
                          power_K3 = mean(raw$p_K3[sel] <= alpha),
                          holm_family1 = mean(raw$holm_family1[sel]),
                          holm_other = mean(raw$holm_other[sel]),
                          power_joint = mean(raw$p_joint_K2to5[sel] <= alpha))
}
summary_table <- do.call(rbind, rows)
write.csv(summary_table, "results/sim_10_semisynthetic_summary.csv", row.names = FALSE)
print(summary_table, digits = 3)
