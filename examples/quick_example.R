################################################################################
#  quick_example.R  --  a minimal use of FCST (runs in well under a minute)
#
#  Run from the folder examples/:   source("quick_example.R")
#  Output: printed results and a figure (quick_example.pdf, also shown on
#  screen in an interactive session).
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
res <- fcst_test(X, labels, scales = radii, B = 999, keep_null = TRUE)
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

## ---- 6. figure ------------------------------------------------------------------
## (a) true types and the hotspot of type 1; (b) estimated types;
## (c) standardised same-type statistics Z_k(r) by radius, for each estimated
##     type, with the Holm-adjusted p-values; (d) permutation distribution of the
##     global statistic G and its observed value.
type_cols <- c("#D55E00", "#0072B2", "#009E73")
type_pch  <- c(16, 17, 15)

pairs <- build_pair_list(X, radii, kernel = "ball")
Z_obs <- standardise_pair(pair_statistic(labels, pairs, K),
                          all_null_moments(labels, pairs, K), K)   # K x radii

draw_quick_figure <- function() {
  op <- par(mfrow = c(2, 2), mar = c(3.2, 3.4, 2, 0.6), mgp = c(2, 0.6, 0),
            cex = 0.9, las = 1)
  g <- seq(0, 1, length.out = 101)
  h <- matrix(hotspot_surface(as.matrix(expand.grid(g, g)), "one_bump"), 101)

  # (a) true types
  plot(X, asp = 1, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n",
       main = "(a) True types")
  contour(g, g, h, levels = c(0.25, 0.75), add = TRUE, drawlabels = FALSE,
          lty = 2, col = "grey45")
  points(X, col = type_cols[true_types], pch = type_pch[true_types], cex = 0.7)

  # (b) estimated types (cluster numbers are arbitrary: the colours of the
  #     estimated types need not match those of the true types)
  plot(X, asp = 1, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n",
       main = sprintf("(b) Estimated types (ARI %.2f)",
                      adjusted_rand_index(true_types, labels)))
  points(X, col = type_cols[labels], pch = type_pch[labels], cex = 0.7)

  # (c) Z_k(r) by radius
  plot(NA, xlim = range(radii), ylim = range(c(0, Z_obs), na.rm = TRUE),
       xlab = "radius r", ylab = expression(Z[k](r)), main = "(c) Same-type pairs")
  abline(h = 0, col = "grey60")
  for (k in 1:K) {
    lines(radii, Z_obs[k, ], col = type_cols[k], lwd = 2)
    points(radii, Z_obs[k, ], col = type_cols[k], pch = type_pch[k])
  }
  legend("topleft", bty = "n", cex = 0.85,
         legend = sprintf("type %d: Holm p = %.3f", 1:K, res$table$p_holm),
         col = type_cols, pch = type_pch, lwd = 2)

  # (d) permutation distribution of G
  hist(res$G_null[, 1], breaks = 40, col = "grey85", border = "white",
       xlim = range(c(res$G_null[, 1], res$G_obs)), xlab = "G under permutation",
       main = sprintf("(d) Global test (p = %.3f)", res$p_global))
  abline(v = res$G_obs, col = "black", lwd = 2.5)
  legend("topright", bty = "n", cex = 0.85, legend = "observed G",
         col = "black", lwd = 2.5)
  par(op)
}

pdf("quick_example.pdf", width = 7, height = 6.4)
draw_quick_figure()
dev.off()
if (interactive()) draw_quick_figure()
cat("\nfigure written to quick_example.pdf\n")
