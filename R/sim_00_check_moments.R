################################################################################
#  sim_00_check_moments.R
#
#  Sanity check of Lemma S1 (Supplementary Material): the closed-form permutation mean and
#  variance of U_k(s) must match the Monte Carlo mean and variance obtained
#  by permuting the labels many times. Checked for both kernels (ball and
#  Gaussian) and both weighting schemes, since the proposition holds for any
#  symmetric pair weights.
#
#  Run from the folder containing fcst_methods.R:
#     source("sim_00_check_moments.R")
################################################################################

source("fcst_methods.R")
set.seed(20261029)

## ---- settings ----------------------------------------------------------------
n       <- 150                        # number of points
K       <- 3                          # number of types
n_perm  <- 20000                      # permutations for the Monte Carlo check
scales_by_kernel <- list(ball = c(0.05, 0.10, 0.20),
                         gaussian = c(1, 2, 4) / 32)

## ---- one pattern and one labelling ------------------------------------------
X      <- sim_points(n, "ipp")
labels <- sample(1:K, n, replace = TRUE)

for (kernel in names(scales_by_kernel)) {
  scales <- scales_by_kernel[[kernel]]
  for (w in c("none", "inverse_intensity")) {
    pairs <- build_pair_list(X, scales, weights = w, kernel = kernel)
    U_perm <- array(NA_real_, dim = c(n_perm, K, length(scales)))
    for (b in 1:n_perm) {
      U_perm[b, , ] <- pair_statistic(labels[sample.int(n)], pairs, K)
    }
    cat("\n=== kernel:", kernel, " weights:", w, "===\n")
    for (k in 1:K) {
      n_k <- sum(labels == k)
      mom <- pair_null_moments(n, n_k, pairs)
      mc_mean <- colMeans(U_perm[, k, ])
      mc_var  <- apply(U_perm[, k, ], 2, var)
      out <- data.frame(scale = scales,
                        mean_exact = mom$mean, mean_MC = mc_mean,
                        var_exact = mom$var, var_MC = mc_var,
                        ratio_var = mc_var / mom$var)
      cat("type", k, "(n_k =", n_k, ")\n")
      print(round(out, 4))
    }
  }
}
## Expected: mean_exact ~ mean_MC and ratio_var ~ 1 (up to Monte Carlo error,
## about +/- 2 * sqrt(2 / n_perm) = +/- 0.02 for the variance ratio).
