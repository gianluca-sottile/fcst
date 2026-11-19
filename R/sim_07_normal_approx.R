################################################################################
#  sim_07_normal_approx.R
#
#  Accuracy of the normal approximation based on the closed-form moments of
#  Lemma S1 (Section S4 of the Supplementary Material): for a single type and a single
#  radius, p_norm = 1 - Phi(Z_k(r)) is compared with the permutation p-value.
#  Under H0 we record the rejection rate of the normal approximation at
#  alpha = 0.05 and 0.01, for several n, radii and type proportions.
#
#  Output: results/sim_07_normal_approx_summary.csv
################################################################################

source("fcst_methods.R")
dir.create("results", showWarnings = FALSE)

## ---- settings ----------------------------------------------------------------
n_rep    <- 2000             # paper: 2000 (cheap: no permutations needed)
n_cores  <- 4
sizes    <- c(100, 300, 1000)
radii    <- c(0.05, 0.10, 0.20)
props    <- c(0.1, 0.33)     # proportion of type-1 points
seed_base <- 7e7

run_one <- function(n, prop) {
  X <- sim_points(n, "hpp")
  labels <- 1 + (runif(n) > prop)             # type 1 with probability prop
  pairs <- build_pair_list(X, radii)
  U <- pair_statistic(labels, pairs, 2)
  mom <- pair_null_moments(n, sum(labels == 1), pairs)
  z <- (U[1, ] - mom$mean) / sqrt(mom$var)
  row <- data.frame(n = n, prop = prop)
  for (m in seq_along(radii)) {
    row[[paste0("p_r", m)]] <- 1 - pnorm(z[m])
  }
  return(row)
}

all_rows <- list()
setting_id <- 0
for (n in sizes) {
  for (prop in props) {
    setting_id <- setting_id + 1
    cat("n =", n, " prop =", prop, "\n")
    rows <- run_replicates(function() run_one(n, prop), n_rep, setting_id,
                           n_cores, seed_base)
    all_rows <- c(all_rows, rows)
  }
}
raw <- do.call(rbind, all_rows)

summary_rows <- list()
counter <- 0
for (n in sizes) {
  for (prop in props) {
    sel <- raw$n == n & raw$prop == prop
    for (m in seq_along(radii)) {
      counter <- counter + 1
      p <- raw[sel, paste0("p_r", m)]
      summary_rows[[counter]] <- data.frame(n = n, prop = prop,
                                            radius = radii[m],
                                            rej_05 = mean(p <= 0.05),
                                            rej_01 = mean(p <= 0.01))
    }
  }
}
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, "results/sim_07_normal_approx_summary.csv",
          row.names = FALSE)
print(summary_table, digits = 3)
## Expected: close to 0.05 / 0.01 for large n and moderate radii; the
## approximation is anti-conservative for small n_k and small radii
## (skewed distribution of U), which is why permutation calibration is the
## default.
