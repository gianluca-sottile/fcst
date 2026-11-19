################################################################################
#  fcst_competitors.R
#
#  Competing tests used in the simulation study (Section 6.6 and Supplementary Section S5 of the paper).
#  Each test is implemented in base R so that all methods share the same data,
#  the same permutation scheme (random labelling: locations fixed, marks
#  permuted) and the same p-value formula (1 + #{null >= obs}) / (B + 1).
#  This makes the comparison fair and avoids differences due to package
#  defaults; differences with the reference implementations are noted below.
#
#  Tests on the TYPE LABELS (need a typing, like FCST)
#    segregation_test()   Diggle, Zheng & Durr (2005): kernel estimates of the
#                         type probabilities and Monte Carlo test of their
#                         constancy (as spatstat::segregation.test).
#    multinomial_scan()   Jung, Kulldorff & Richard (2010): circular scan
#                         windows, multinomial likelihood ratio.
#
#  Tests on the CURVES (no typing)
#    functional_scan()    circular scan windows, L2 distance between the mean
#                         representation inside and outside the window, in the
#                         spirit of Frevent et al. (2021).
#    mark_weighted_K()    global mark-weighted K-function with functional
#                         marks (Comas et al. 2011; Ghorbani et al. 2021;
#                         global counterpart of D'Angelo et al. 2024).
#    dcov_test()          distance covariance between locations and curves
#                         (Szekely-Rizzo; Lyons 2013).
#
#  All functions return list(statistic = ..., p_value = ...).
#  Requires fcst_methods.R to be sourced first (build_pair_list,
#  draw_permutation, fpca_scores, default_bandwidth).
################################################################################


## =============================================================================
## Helpers
## =============================================================================

perm_p_value <- function(obs, null) {
  # Upper-tail Monte Carlo p-value, valid with randomly drawn permutations.
  return((1 + sum(null >= obs)) / (length(null) + 1))
}

row_cumsum <- function(M) {
  # Cumulative sums along each row of a matrix (column by column, so that
  # the loop is over the m columns and each step is a vector operation).
  if (ncol(M) > 1) {
    for (j in 2:ncol(M)) {
      M[, j] <- M[, j] + M[, j - 1]
    }
  }
  return(M)
}

neighbour_order <- function(X, max_frac = 0.5) {
  # For every point (window centre) the indices of the points ordered by
  # distance, truncated at max_frac * n: row i, column c = the c-th nearest
  # point to point i (the point itself is in column 1). Circular scan windows
  # are the sets {first c points of row i}.
  n <- nrow(X)
  m <- floor(max_frac * n)
  D <- as.matrix(dist(X))
  O <- matrix(0L, nrow = n, ncol = m)
  for (i in 1:n) {
    O[i, ] <- order(D[i, ])[1:m]
  }
  return(O)
}

xlogx_ratio <- function(x, y) {
  # x * log(x / y) with the convention 0 * log(0) = 0 (vectorised).
  out <- x * log(x / y)
  out[x == 0] <- 0
  return(out)
}


## =============================================================================
## 1. Segregation test (Diggle, Zheng & Durr 2005)
## =============================================================================
## p_k(x_i) is estimated by leave-one-out Nadaraya-Watson smoothing of the type
## indicators with a Gaussian kernel; the statistic is
##   T = sum_i sum_k { p_hat_k(x_i) - n_k / n }^2,
## large when type probabilities vary in space. spatstat::segregation.test
## uses the same statistic; it selects the bandwidth by cross-validation,
## whereas here the bandwidth is fixed (normal-reference rule) for speed.

segregation_test <- function(X, labels, B = 999, bw = NULL, strata = NULL,
                             bw_grid = NULL,
    progress = FALSE, progress_label = "Segregation") {
  # bw: NULL      normal-reference rule, fixed;
  #     a number  fixed bandwidth;
  #     "cv"      likelihood cross-validation on the OBSERVED labels, then
  #               kept fixed across permutations (as spatstat's
  #               segregation.test does by default);
  #     "cv_perm" cross-validation repeated for every permuted labelling,
  #               so that the selection step is part of the statistic and
  #               the test is exact (slower).
  n <- nrow(X)
  D2 <- as.matrix(dist(X))^2
  cats <- sort(unique(labels))
  p_bar <- as.numeric(table(factor(labels, levels = cats))) / n
  stat_given_W <- function(lab, Wk) {
    row_tot <- rowSums(Wk)
    total <- 0
    for (j in seq_along(cats)) {
      p_hat <- as.vector(Wk %*% as.numeric(lab == cats[j])) / row_tot
      total <- total + sum((p_hat - p_bar[j])^2)
    }
    return(total)
  }
  kernel_matrix <- function(h) {
    Wk <- exp(-D2 / (2 * h^2))
    diag(Wk) <- 0
    return(Wk)
  }
  if (identical(bw, "cv_perm")) {
    if (is.null(bw_grid)) {
      bw_grid <- default_bandwidth(X) * exp(seq(log(0.25), log(4), length.out = 15))
    }
    W_list <- list()
    for (g in seq_along(bw_grid)) W_list[[g]] <- kernel_matrix(bw_grid[g])
    statistic <- function(lab) {
      same <- outer(lab, lab, "==")
      best <- 1
      best_val <- -Inf
      for (g in seq_along(W_list)) {
        p_own <- rowSums(W_list[[g]] * same) / rowSums(W_list[[g]])
        val <- sum(log(pmax(p_own, 1e-10)))
        if (val > best_val) {
          best_val <- val
          best <- g
        }
      }
      return(stat_given_W(lab, W_list[[best]]))
    }
  } else {
    if (is.null(bw)) bw <- default_bandwidth(X)
    if (identical(bw, "cv")) bw <- cv_bandwidth_types(X, labels, bw_grid)
    Wk <- kernel_matrix(bw)
    statistic <- function(lab) stat_given_W(lab, Wk)
  }
  obs <- statistic(labels)
  null <- numeric(B)
  t0 <- Sys.time()
  for (b in 1:B) {
    if (progress) progress_report(b, B, t0, progress_label)
    null[b] <- statistic(labels[draw_permutation(n, strata)])
  }
  return(list(statistic = obs, p_value = perm_p_value(obs, null), bw = bw))
}


## =============================================================================
## 2. Multinomial spatial scan statistic (Jung, Kulldorff & Richard 2010)
## =============================================================================
## Windows: the c nearest points to each point, c = 2, ..., floor(max_frac n).
## For a window with counts C_k (size c) and totals N_k:
##   LLR = sum_k [ C_k log(C_k / c) + (N_k - C_k) log((N_k - C_k) / (n - c)) ]
##         - sum_k N_k log(N_k / n).
## The scan statistic is the maximum LLR over windows; it is calibrated by
## random labelling, as in SaTScan.

multinomial_scan <- function(X, labels, B = 999, max_frac = 0.5,
                             O = NULL, strata = NULL,
    progress = FALSE, progress_label = "Multinomial scan") {
  n <- nrow(X)
  if (is.null(O)) O <- neighbour_order(X, max_frac)
  m <- ncol(O)
  sizes <- matrix(rep(1:m, each = n), nrow = n)       # window size c
  cats <- sort(unique(labels))
  N_k <- as.numeric(table(factor(labels, levels = cats)))
  constant <- sum(xlogx_ratio(N_k, n))

  scan_stat <- function(lab) {
    L <- matrix(lab[O], nrow = n)                      # labels by window
    llr <- matrix(-constant, nrow = n, ncol = m)
    for (j in seq_along(cats)) {
      C <- row_cumsum(1 * (L == cats[j]))
      llr <- llr + xlogx_ratio(C, sizes) +
        xlogx_ratio(N_k[j] - C, n - sizes)
    }
    llr[, 1] <- -Inf                                   # windows of size >= 2
    best <- which.max(llr)
    return(list(value = llr[best],
                centre = ((best - 1) %% n) + 1,
                size = ((best - 1) %/% n) + 1))
  }

  obs <- scan_stat(labels)
  null <- numeric(B)
  t0 <- Sys.time()
  for (b in 1:B) {
    if (progress) progress_report(b, B, t0, progress_label)
    null[b] <- scan_stat(labels[draw_permutation(n, strata)])$value
  }
  return(list(statistic = obs$value, p_value = perm_p_value(obs$value, null),
              window_centre = obs$centre, window_size = obs$size))
}


## =============================================================================
## 3. Functional scan statistic (curves, no typing)
## =============================================================================
## Windows as in multinomial_scan(). Curves are represented by their FPCA
## scores S (n x q, centred). For a window W of size c,
##   F(W) = n / (c (n - c)) * || sum_{i in W} s_i ||^2
##        = c (n - c) / n * || mean_W - mean_outside ||^2,
## the between-group sum of squares of a one-way functional ANOVA with
## groups "inside" and "outside". The scan statistic is max_W F(W),
## calibrated by permuting the curves over the locations. This follows the
## idea of the parametric functional scan of Frevent et al. (2021); it is our
## implementation, not their software (R package HDSpatialScan).

functional_scan <- function(X, Z, B = 999, max_frac = 0.5,
                            var_explained = 0.90, O = NULL, strata = NULL,
    progress = FALSE, progress_label = "Functional scan") {
  n <- nrow(X)
  if (is.null(O)) O <- neighbour_order(X, max_frac)
  m <- ncol(O)
  S <- fpca_scores(Z, var_explained)                  # centred scores
  sizes <- matrix(rep(1:m, each = n), nrow = n)
  scale_c <- n / (sizes * (n - sizes))

  scan_stat <- function(S_perm) {
    sq_norm <- matrix(0, nrow = n, ncol = m)
    for (q in 1:ncol(S_perm)) {
      cum <- row_cumsum(matrix(S_perm[O, q], nrow = n))
      sq_norm <- sq_norm + cum^2
    }
    F_stat <- scale_c * sq_norm
    F_stat[, 1] <- -Inf
    return(max(F_stat))
  }

  obs <- scan_stat(S)
  null <- numeric(B)
  t0 <- Sys.time()
  for (b in 1:B) {
    if (progress) progress_report(b, B, t0, progress_label)
    null[b] <- scan_stat(S[draw_permutation(n, strata), , drop = FALSE])
  }
  return(list(statistic = obs, p_value = perm_p_value(obs, null)))
}


## =============================================================================
## 4. Global mark-weighted K-function with functional marks
## =============================================================================
## K_f(r) = sum_{i != j} 1{d_ij <= r} w_ij * || z_i - z_j ||^2
## (squared L2 distance between representations as the mark function).
## Small values mean that nearby points carry similar curves. Each radius is
## standardised by the mean and s.d. of the B + 1 values (observed included,
## which preserves exchangeability), and the statistic is
##   T = max_r { -(K_f(r) - mean) / sd },
## a max-type global test across radii in the spirit of global envelope
## tests (Myllymaki et al. 2017).

mark_weighted_K <- function(X, Z, radii, B = 999, weights = "none",
                            pairs = NULL, strata = NULL,
    progress = FALSE, progress_label = "Mark-weighted K") {
  n <- nrow(X)
  if (is.null(pairs)) {
    pairs <- build_pair_list(X, radii, weights = weights, kernel = "ball")
  }
  R <- length(radii)
  G <- tcrossprod(Z)
  sq <- diag(G)
  Dmark <- outer(sq, sq, "+") - 2 * G                   # squared L2 distances

  K_curve <- function(perm) {
    d_pair <- Dmark[cbind(perm[pairs$i], perm[pairs$j])]
    return(2 * colSums(pairs$a * d_pair))
  }

  all_K <- matrix(NA_real_, nrow = B + 1, ncol = R)
  all_K[1, ] <- K_curve(1:n)
  t0 <- Sys.time()
  for (b in 1:B) {
    if (progress) progress_report(b, B, t0, progress_label)
    all_K[b + 1, ] <- K_curve(draw_permutation(n, strata))
  }
  mu <- colMeans(all_K)
  sdv <- apply(all_K, 2, sd)
  Zs <- -(all_K - matrix(mu, nrow = B + 1, ncol = R, byrow = TRUE)) /
    matrix(sdv, nrow = B + 1, ncol = R, byrow = TRUE)
  T_all <- apply(Zs, 1, max)
  return(list(statistic = T_all[1], p_value = perm_p_value(T_all[1], T_all[-1])))
}


## =============================================================================
## 5. Distance covariance between locations and curves
## =============================================================================
## V-statistic dCov^2 = mean(A * Bc), where A and Bc are the double-centred
## Euclidean distance matrices of the locations and of the representations.
## Consistent against every form of dependence between curves and locations
## (Lyons 2013), but it does not say which shapes are involved.

double_centre <- function(D) {
  rm <- rowMeans(D)
  cm <- colMeans(D)
  return(D - outer(rm, cm, "+") + mean(D))
}

dcov_test <- function(X, Z, B = 999, strata = NULL,
    progress = FALSE, progress_label = "Distance covariance") {
  n <- nrow(X)
  A <- double_centre(as.matrix(dist(X)))
  Bc <- double_centre(as.matrix(dist(Z)))
  obs <- mean(A * Bc)
  null <- numeric(B)
  t0 <- Sys.time()
  for (b in 1:B) {
    if (progress) progress_report(b, B, t0, progress_label)
    perm <- draw_permutation(n, strata)
    null[b] <- mean(A * Bc[perm, perm])
  }
  return(list(statistic = obs, p_value = perm_p_value(obs, null)))
}
