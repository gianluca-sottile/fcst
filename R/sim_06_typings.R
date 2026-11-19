################################################################################
#  sim_06_typings.R
#
#  FCST with several typing (clustering) algorithms and joint calibration
#  (Section 3.4 and Proposition 1 of the paper).
#
#  Typings: all algorithms available on this machine among
#    fpca_kmeans  FPCA scores + k-means
#    fpcac        FPCAC: B-spline FPCA + trimmed k-means (clustEff::fpcac)
#    fpca_mclust  FPCA scores + Gaussian mixture (package mclust)
#    fpca_ward    FPCA scores + Ward agglomerative clustering
#    fpca_pam     FPCA scores + partitioning around medoids (package cluster)
#    clusteff     shape-distance clustering (clustEff, Sottile & Adelfio 2019),
#                 distance component only (see check_clusteff.R)
#  All typings use the same Fourier-modulus representation, K = 3, and are
#  computed on the same curves. SNR = 0.05, so that the typings are imperfect
#  and disagree with one another.
#
#  For each replicate we record
#    ari_<alg>       adjusted Rand index of each typing with the true types
#    p_<alg>         p-value of FCST with that typing alone
#    p_naive_min     min over typings of p_<alg>: "pick the clustering that
#                    works" -- NOT a valid p-value
#    p_bonferroni    min(1, A * p_naive_min)
#    p_joint         joint max-T p-value over all typings (Proposition 1)
#
#  Output: results/sim_06_typings_raw.csv, results/sim_06_typings_summary.csv
################################################################################

source("fcst_methods.R")
dir.create("results", showWarnings = FALSE)

## ---- settings ----------------------------------------------------------------
n_rep    <- 500              # paper: 1000 for beta = 0, 500 otherwise
B        <- 999              # paper: 999
n_cores  <- 1
n        <- 300
K        <- 3
snr      <- 0.05
radii    <- c(0.05, 0.10, 0.15, 0.20)
alpha    <- 0.05
patterns <- c("one_bump", "band")
betas    <- c(0, 3, 6)

seed_base <- 6e7
set.seed(seed_base)

## ---- which algorithms can be used? -----------------------------------------------
## Each candidate is tried once on a small example; algorithms that are not
## installed or fail are dropped with a message.
## All typings act on the same FPCA scores of the Fourier modulus (fair
## comparison), except clustEff, which needs the curves themselves.
candidates <- c("fpca_kmeans", "fpca_mclust", "fpca_ward", "fpca_pam",
                "fpcac", "clusteff")
## fpcac = FPCAC of clustEff (B-spline FPCA + trimmed k-means, 5% trimming);
## clusteff = shape-distance clustering, distance component only (see
## check_clusteff.R). Both act on the Fourier modulus.
test_types <- rep(1:K, each = 20)
test_Y <- sim_curves(test_types, snr = 1)
algorithms <- character(0)
for (alg in candidates) {
  ok <- tryCatch({
    lab <- type_curves(test_Y, K, representation = "spectrum", algorithm = alg)
    length(lab) == length(test_types)
  }, error = function(e) {
    cat("algorithm", alg, "dropped:", conditionMessage(e), "\n")
    FALSE
  })
  if (ok) algorithms <- c(algorithms, alg)
}
A <- length(algorithms)
cat("typing algorithms used:", algorithms, "\n")

## ---- one replicate ----------------------------------------------------------
run_one <- function(pattern, beta) {
  X <- sim_points(n, "hpp")
  types <- sim_types(type_probabilities(X, K = K, beta = beta,
                                        pattern = pattern))
  Y <- sim_curves(types, snr = snr)
  label_matrix <- matrix(0L, nrow = n, ncol = A)
  for (a in 1:A) {
    label_matrix[, a] <- type_curves(Y, K, representation = "spectrum",
                                     algorithm = algorithms[a])
  }
  test <- fcst_test(X, label_matrix, radii, B = B)
  row <- data.frame(pattern = pattern, beta = beta)
  for (a in 1:A) {
    row[[paste0("ari_", algorithms[a])]] <-
      adjusted_rand_index(types, label_matrix[, a])
  }
  for (a in 1:A) {
    row[[paste0("p_", algorithms[a])]] <- test$p_typing[a]
  }
  row$p_naive_min  <- min(test$p_typing)
  row$p_bonferroni <- min(1, A * min(test$p_typing))
  row$p_joint      <- test$p_global
  return(row)
}

## ---- loop over settings -------------------------------------------------------
all_rows <- list()
setting_id <- 0
for (pattern in patterns) {
  for (beta in betas) {
    setting_id <- setting_id + 1
    cat("setting", setting_id, ": pattern =", pattern, " beta =", beta, "\n")
    rows <- run_replicates(function() run_one(pattern, beta), n_rep,
                           setting_id, n_cores, seed_base)
    all_rows <- c(all_rows, rows)
  }
}
raw <- do.call(rbind, all_rows)
write.csv(raw, "results/sim_06_typings_raw.csv", row.names = FALSE)

## ---- summary ----------------------------------------------------------------------
settings <- unique(raw[, c("pattern", "beta")])
summary_rows <- list()
for (s in 1:nrow(settings)) {
  sel <- raw$pattern == settings$pattern[s] & raw$beta == settings$beta[s]
  row <- data.frame(settings[s, ], n_rep = sum(sel))
  for (alg in algorithms) {
    row[[paste0("ari_", alg)]] <- mean(raw[sel, paste0("ari_", alg)])
  }
  for (alg in algorithms) {
    row[[paste0("rej_", alg)]] <- mean(raw[sel, paste0("p_", alg)] <= alpha)
  }
  row$rej_naive_min  <- mean(raw$p_naive_min[sel] <= alpha)
  row$rej_bonferroni <- mean(raw$p_bonferroni[sel] <= alpha)
  row$rej_joint      <- mean(raw$p_joint[sel] <= alpha)
  summary_rows[[s]] <- row
}
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, "results/sim_06_typings_summary.csv", row.names = FALSE)
print(summary_table, digits = 3)
## Expected at beta = 0: every single typing and the joint test near alpha;
## rej_naive_min above alpha (selection after the fact); Bonferroni below
## alpha. At beta > 0: the joint test close to the best single typing.
