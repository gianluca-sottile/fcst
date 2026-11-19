################################################################################
#  quick_example.R  --  a minimal use of FCST (runs in well under a minute)
#
#  Run from the folder examples/:   source("quick_example.R")
################################################################################

source(file.path("..", "R", "fcst_methods.R"))
set.seed(2026)

## ---- 1. data: 300 events, three curve types, type 1 concentrated in a hotspot --------
n <- 300
K <- 3
X <- sim_points(n, "hpp")                                  # locations in [0, 1]^2
P <- type_probabilities(X, K = K, beta = 4, pattern = "one_bump")
true_types <- sim_types(P)
Y <- sim_curves(true_types, snr = 0.05)                    # noisy transient curves

## ---- 2. clustering of the curves (the locations are NOT used) -----------------------
labels <- type_curves(Y, K, representation = "spectrum", algorithm = "fpca_kmeans")
cat("adjusted Rand index with the true types:",
    round(adjusted_rand_index(true_types, labels), 2), "\n")

## ---- 3. the test -----------------------------------------------------------------
radii <- c(0.05, 0.10, 0.15, 0.20)
res <- fcst_test(X, labels, scales = radii, B = 999)
cat("\nglobal test of random labelling: p =", res$p_global, "\n")
cat("type-specific tests (p_holm: Holm-adjusted p-values):\n")
print(res$table[, c("type", "n_k", "T_obs", "p_raw", "p_holm")], digits = 3)

## ---- 4. two clusterings calibrated jointly ------------------------------------------
if (requireNamespace("mclust", quietly = TRUE)) {
  labels_mix <- type_curves(Y, K, representation = "spectrum", algorithm = "fpca_mclust")
  both <- cbind(kmeans = labels, mixture = labels_mix)
  res2 <- fcst_test(X, both, scales = radii, B = 999)
  cat("\njoint test over the two clusterings: p =", res2$p_global, "\n")
  cat("each clustering separately:", round(res2$p_typing, 4), "\n")
}

## ---- 5. covariate-adjusted test: permutations within strata of a covariate ------------
z <- X[, 2]                                               # here, the second coordinate
res3 <- fcst_test(X, labels, scales = radii, B = 999, strata = make_strata(z, n_bins = 4))
cat("\ntest restricted to 4 strata of z: p =", res3$p_global, "\n")
