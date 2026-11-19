################################################################################
#  fcst_methods.R
#
#  Functional Clustering Spatial Test (FCST): implementation of the methodology
#  described in the manuscript the accompanying paper.
#
#  Authors : Gianluca Sottile, Marco D'Antoni, Giada Adelfio
#  Version : 1.0
#
#  Design principles
#  -----------------
#  * Base R + 'stats' only for the core. Optional packages are used only if
#    installed ('cluster' for the silhouette, 'tclust' for trimmed k-means).
#  * Plain syntax: explicit for-loops and small functions, no pipes, no
#    non-standard evaluation, so that every step can be inspected with
#    browser() or print().
#  * Every function returns a plain list / vector / matrix.
#
#  Contents
#  --------
#   1. Simulation of point patterns          (sim_points_*)
#   2. Simulation of types and curves        (type_probabilities, sim_types,
#                                             sim_curves)
#   3. Curve representations                 (represent_curves)
#   4. Functional clustering ("typing")      (fpca_scores, type_curves,
#                                             relabel_canonical)
#   5. Spatial ingredients                   (kde_intensity, kde_at,
#                                             build_pair_list)
#   6. Test statistics                       (pair_statistic,
#                                             pair_null_moments,
#                                             centroid_statistic)
#   7. Permutation machinery                 (draw_permutation)
#   8. The test                              (fcst_test)
#   9. Misclassification helpers             (add_label_noise,
#                                             attenuation_gamma, delta_signal)
#  10. Local (point-wise) diagnostics        (local_same_type)
#
#  Notation used in the comments (same as in the paper)
#  ----------------------------------------------------
#   n        number of points
#   X        n x 2 matrix of coordinates
#   Y        n x T matrix of curves evaluated on a common grid
#   c        vector of n type labels (0 = unassigned / trimmed, 1..K types)
#   s_1 < ... < s_S  spatial scales (Gaussian bandwidths, or radii for the
#            ball kernel)
#   a_ij(s)  = kappa_s(d_ij) * w_ij, kappa_s(d) = exp(-d^2 / (2 s^2)) or
#              1{d <= s};  w_ij = 1 or inverse-intensity weights
#   U_k(s)   = sum_{i != j} 1{c_i = k} 1{c_j = k} a_ij(s)     (ordered pairs)
################################################################################

# Optional packages (tclust, mclust, cluster, clustEff) are called with "::"
# inside the functions that need them; the core of the method uses base R only.


## =============================================================================
## 1. SIMULATION OF POINT PATTERNS (unit square window, n fixed)
## =============================================================================
## All generators work conditionally on the number of points n, which is the
## natural conditioning for a random-labelling test.

sim_points_hpp <- function(n) {
  # Homogeneous Poisson process conditioned on N = n: iid uniform points.
  X <- cbind(runif(n), runif(n))
  colnames(X) <- c("x", "y")
  return(X)
}

sim_points_ipp <- function(n, intensity_fun) {
  # Inhomogeneous Poisson process conditioned on N = n: iid points with
  # density proportional to intensity_fun(x, y). Rejection sampling.
  # intensity_fun must be vectorised and bounded on the unit square.
  grid_check <- expand.grid(x = seq(0, 1, length.out = 101),
                            y = seq(0, 1, length.out = 101))
  lambda_max <- max(intensity_fun(grid_check$x, grid_check$y)) * 1.05
  X <- matrix(NA_real_, nrow = 0, ncol = 2)
  while (nrow(X) < n) {
    m <- 2 * (n - nrow(X))
    cand <- cbind(runif(m), runif(m))
    keep <- runif(m) < intensity_fun(cand[, 1], cand[, 2]) / lambda_max
    X <- rbind(X, cand[keep, , drop = FALSE])
  }
  X <- X[1:n, , drop = FALSE]
  colnames(X) <- c("x", "y")
  return(X)
}

sim_points_thomas <- function(n, n_parents = 15, sigma = 0.03,
                              fixed_parents = NULL) {
  # Thomas-type cluster process conditioned on N = n (approximation):
  # parents uniform on the square; each offspring picks a parent at random
  # and is displaced by an isotropic Gaussian; offspring outside the window
  # are discarded and redrawn.
  # fixed_parents (optional, m x 2 matrix): parents placed at given
  # locations (e.g. at the hotspot centres, so that clusters of events and
  # concentration of a type overlap); the other n_parents - m are uniform.
  n_fixed <- 0
  if (!is.null(fixed_parents)) n_fixed <- nrow(fixed_parents)
  n_free <- max(n_parents - n_fixed, 0)
  parents <- cbind(runif(n_free), runif(n_free))
  if (n_fixed > 0) parents <- rbind(as.matrix(fixed_parents), parents)
  n_par <- nrow(parents)
  X <- matrix(NA_real_, nrow = 0, ncol = 2)
  while (nrow(X) < n) {
    m <- 2 * (n - nrow(X))
    which_parent <- sample.int(n_par, m, replace = TRUE)
    cand <- parents[which_parent, , drop = FALSE] +
      matrix(rnorm(2 * m, sd = sigma), ncol = 2)
    inside <- cand[, 1] >= 0 & cand[, 1] <= 1 & cand[, 2] >= 0 & cand[, 2] <= 1
    X <- rbind(X, cand[inside, , drop = FALSE])
  }
  X <- X[1:n, , drop = FALSE]
  colnames(X) <- c("x", "y")
  return(X)
}

sim_points <- function(n, process = "hpp", ...) {
  # Wrapper: process in {"hpp", "ipp", "thomas"}.
  # For "ipp" the default intensity is the linear trend lambda(x, y) ~ x.
  if (process == "hpp") {
    return(sim_points_hpp(n))
  }
  if (process == "ipp") {
    args <- list(...)
    if (is.null(args$intensity_fun)) {
      f <- function(x, y) 0.1 + x
    } else {
      f <- args$intensity_fun
    }
    return(sim_points_ipp(n, f))
  }
  if (process == "thomas") {
    return(sim_points_thomas(n, ...))
  }
  stop("Unknown process: ", process)
}


## =============================================================================
## 2. SIMULATION OF TYPES AND FUNCTIONAL MARKS
## =============================================================================

bump_profile <- function(X, centre, scale, profile = "gauss", df = 3) {
  # Radial bump centred at 'centre'. Both profiles equal 2 pi s^2 times a
  # bivariate density with mean 'centre' and covariance s^2 I, so they have
  # the same integral (2 pi s^2) and the same spread (second moment):
  #   "gauss": exp(-d^2 / (2 s^2))                                  (peak 1)
  #   "t"    : df / (df - 2) * (1 + d^2 / ((df - 2) s^2))^(-(df + 2) / 2)
  #            bivariate Student t with df degrees of freedom   (peak df/(df-2))
  # With the same total excess and the same spread, the t profile has a
  # denser core and heavier tails: "denser hotspots". Requires df > 2.
  d2 <- (X[, 1] - centre[1])^2 + (X[, 2] - centre[2])^2
  if (profile == "gauss") {
    return(exp(-d2 / (2 * scale^2)))
  }
  if (profile == "t") {
    return(df / (df - 2) * (1 + d2 / ((df - 2) * scale^2))^(-(df + 2) / 2))
  }
  stop("Unknown profile: ", profile)
}

hotspot_centres <- function(pattern) {
  # Centres of the hotspots of type 1 (NULL for "none" and "band").
  base <- sub("_t$", "", pattern)
  if (base == "one_bump") return(matrix(c(0.3, 0.3), ncol = 2))
  if (base == "two_bumps") return(matrix(c(0.25, 0.25, 0.75, 0.75),
                                         ncol = 2, byrow = TRUE))
  return(NULL)
}

hotspot_surface <- function(X, pattern = "one_bump") {
  # Spatial "excess" surface h(x) >= 0 used for type 1 (in [0, 1] except
  # for the t profiles, whose peak is df / (df - 2) = 3).
  #   "none"        : h = 0 (null hypothesis)
  #   "one_bump"    : one isotropic Gaussian bump centred at (0.3, 0.3)
  #   "two_bumps"   : two bumps at (0.25, 0.25) and (0.75, 0.75)
  #   "band"        : an elongated band along the line y = x (fault-like)
  #   "one_bump_t", "two_bumps_t": as above with a bivariate Student-t
  #                   profile (df = 3) with the same integral and spread:
  #                   denser hotspots, heavier tails
  x <- X[, 1]
  y <- X[, 2]
  profile <- "gauss"
  if (grepl("_t$", pattern)) profile <- "t"
  base <- sub("_t$", "", pattern)
  if (base == "none") {
    h <- rep(0, nrow(X))
  } else if (base == "one_bump") {
    h <- bump_profile(X, c(0.3, 0.3), 0.10, profile)
  } else if (base == "two_bumps") {
    h1 <- bump_profile(X, c(0.25, 0.25), 0.08, profile)
    h2 <- bump_profile(X, c(0.75, 0.75), 0.08, profile)
    h <- pmax(h1, h2)
  } else if (base == "band") {
    dist_to_line <- abs(y - x) / sqrt(2)
    h <- exp(-dist_to_line^2 / (2 * 0.05^2))
  } else {
    stop("Unknown pattern: ", pattern)
  }
  return(h)
}

type_probabilities <- function(X, K = 3, beta = 0, pattern = "one_bump") {
  # Multinomial-logit type model:
  #   p_k(x) proportional to exp(alpha_k + beta * h_k(x)),
  # with h_1 = hotspot_surface(X, pattern) and h_k = 0 for k >= 2.
  # beta = 0 gives constant type probabilities (null hypothesis).
  h <- hotspot_surface(X, pattern)
  eta <- matrix(0, nrow = nrow(X), ncol = K)
  eta[, 1] <- beta * h
  P <- exp(eta)
  P <- P / rowSums(P)
  return(P)
}

sim_types <- function(P) {
  # Draw one label per row of the n x K probability matrix P.
  n <- nrow(P)
  K <- ncol(P)
  types <- integer(n)
  for (i in 1:n) {
    types[i] <- sample.int(K, 1, prob = P[i, ])
  }
  return(types)
}

type_shape_parameters <- function(K) {
  # Default "seismogram-like" parameters for each type: dominant frequency
  # (cycles on [0, 1]) and time of the envelope peak.
  freq <- seq(6, 6 + 4 * (K - 1), by = 4)          # 6, 10, 14, ...
  peak <- seq(0.15, 0.15 + 0.10 * (K - 1), by = 0.10)  # 0.15, 0.25, ...
  return(list(freq = freq, peak = peak))
}

sim_curves <- function(types, n_grid = 200, snr = 1, amp_sdlog = 0.5,
                       amp_multiplier = NULL, shape = NULL) {
  # Transient, seismogram-like curves on t in [0, 1]:
  #   s_i(t) = A_i * e_k(t) * sin(2 * pi * f_k * t + phi_i)
  #   e_k(t) = (t / tau_k) * exp(1 - t / tau_k)       (peak 1 at tau_k)
  #   A_i    = amp_multiplier_i * exp(N(0, amp_sdlog^2))
  #   phi_i  ~ U(0, 2 pi)                               (random phase)
  #   Y_i(t) = s_i(t) + eps_i(t),  eps_i(t) iid N(0, sigma_i^2)
  #   sigma_i^2 = var_t(s_i) / snr   (SNR = signal variance / noise variance)
  # amp_multiplier (length n) can depend on location: used to create
  # confounding in sim_04. snr can be a scalar or a vector of length n
  # (curve quality depending on location, sim_04b).
  n <- length(types)
  K <- max(types)
  if (is.null(shape)) shape <- type_shape_parameters(K)
  if (is.null(amp_multiplier)) amp_multiplier <- rep(1, n)
  tt <- seq(0, 1, length.out = n_grid)
  Y <- matrix(0, nrow = n, ncol = n_grid)
  for (i in 1:n) {
    k <- types[i]
    tau <- shape$peak[k]
    envelope <- (tt / tau) * exp(1 - tt / tau)
    phase <- runif(1, 0, 2 * pi)
    amplitude <- amp_multiplier[i] * exp(rnorm(1, 0, amp_sdlog))
    signal <- amplitude * envelope * sin(2 * pi * shape$freq[k] * tt + phase)
    snr_i <- snr[1]
    if (length(snr) > 1) snr_i <- snr[i]      # SNR may vary by curve (sim_04b)
    sigma <- sqrt(var(signal) / snr_i)
    Y[i, ] <- signal + rnorm(n_grid, 0, sigma)
  }
  attr(Y, "grid") <- tt
  return(Y)
}


## =============================================================================
## 3. CURVE REPRESENTATIONS
## =============================================================================

represent_curves <- function(Y, method = "spectrum", n_freq = NULL) {
  # Map each curve to a finite-dimensional representation.
  #   "raw"        : the curve itself (sensitive to amplitude and phase)
  #   "unit_norm"  : curve divided by its L2 norm (amplitude-invariant)
  #   "spectrum"   : modulus of the discrete Fourier transform of the
  #                  centred curve, frequencies 1..n_freq, normalised to unit
  #                  norm (amplitude- and phase-invariant)
  n <- nrow(Y)
  if (method == "raw") {
    return(Y)
  }
  if (method == "unit_norm") {
    Z <- Y
    for (i in 1:n) {
      Z[i, ] <- Y[i, ] / sqrt(sum(Y[i, ]^2))
    }
    return(Z)
  }
  if (method == "spectrum") {
    n_grid <- ncol(Y)
    if (is.null(n_freq)) n_freq <- floor(n_grid / 4)
    Z <- matrix(0, nrow = n, ncol = n_freq)
    for (i in 1:n) {
      centred <- Y[i, ] - mean(Y[i, ])
      modulus <- Mod(fft(centred))[2:(n_freq + 1)]
      Z[i, ] <- modulus / sqrt(sum(modulus^2))
    }
    return(Z)
  }
  stop("Unknown representation: ", method)
}


## =============================================================================
## 4. FUNCTIONAL CLUSTERING ("TYPING")
## =============================================================================
## Requirement for validity (Assumption 1 and Proposition 1 in the paper): the typing map must
## use the curves ONLY (never the locations) and must not depend on the order
## of the curves except through randomness independent of the data.

fpca_scores <- function(Z, var_explained = 0.90, max_comp = 20) {
  # Discretised FPCA via SVD of the centred representation matrix.
  # On an equispaced grid this equals L2-FPCA up to a constant factor.
  pca <- prcomp(Z, center = TRUE, scale. = FALSE)
  share <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
  n_comp <- which(share >= var_explained)[1]
  n_comp <- max(2, min(n_comp, max_comp))
  scores <- pca$x[, 1:n_comp, drop = FALSE]
  return(scores)
}

relabel_canonical <- function(labels, Z) {
  # Canonical labelling, so that "type k" is defined by the curves and not by
  # an arbitrary label returned by the algorithm:
  #   types are ordered by decreasing size; ties are broken by the L2 norm
  #   of the mean representation. Label 0 (trimmed) is left untouched.
  ids <- sort(unique(labels[labels > 0]))
  sizes <- numeric(length(ids))
  norms <- numeric(length(ids))
  for (j in seq_along(ids)) {
    members <- labels == ids[j]
    sizes[j] <- sum(members)
    norms[j] <- sqrt(sum(colMeans(Z[members, , drop = FALSE])^2))
  }
  ord <- order(-sizes, -norms)
  new_labels <- labels
  for (j in seq_along(ord)) {
    new_labels[labels == ids[ord[j]]] <- j
  }
  return(new_labels)
}

cluster_kmeans <- function(scores, K, nstart = 25) {
  fit <- kmeans(scores, centers = K, nstart = nstart, iter.max = 100)
  return(fit$cluster)
}

cluster_trimmed_kmeans <- function(scores, K, trim = 0.05) {
  # Trimmed k-means (Garcia-Escudero & Gordaliza 2005) via package 'tclust'.
  # Trimmed curves receive label 0.
  if (!requireNamespace("tclust", quietly = TRUE)) {
    stop("Package 'tclust' is required for trimmed k-means.")
  }
  fit <- tclust::tkmeans(scores, k = K, alpha = trim)
  return(fit$cluster)
}

cluster_ward <- function(Z, K) {
  # Agglomerative clustering with Ward's criterion on L2 distances.
  hc <- hclust(dist(Z), method = "ward.D2")
  return(cutree(hc, k = K))
}

cluster_mclust <- function(scores, K) {
  # Gaussian mixture model on the FPCA scores (package 'mclust'); the
  # covariance structure is chosen by BIC among mclust's default models.
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required for model-based clustering.")
  }
  fit <- mclust::Mclust(scores, G = K, verbose = FALSE)
  return(fit$classification)
}

cluster_pam <- function(Z, K) {
  # Partitioning around medoids on L2 distances (package 'cluster').
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Package 'cluster' is required for PAM.")
  }
  fit <- cluster::pam(dist(Z), k = K, diss = TRUE)
  return(fit$clustering)
}

cluster_clusteff <- function(Z, K, alpha = 0.5, step = "distance") {
  # Shape-distance clustering of Sottile & Adelfio (2019), package clustEff
  # (>= 1.0.0). Z (n x p) is passed as p x n (one column per curve).
  # Settings chosen with check_clusteff.R on the Fourier modulus: the shape
  # component (signs of derivatives) is uninformative on the noise floor of
  # the spectra, so only the distance component is used (step = "distance"):
  # ARI 0.98 at SNR = 1 against 0.66 with the default step = "both".
  if (!requireNamespace("clustEff", quietly = TRUE)) {
    stop("Package 'clustEff' is required.")
  }
  fit <- clustEff::clustEff(t(Z), k = K, alpha = alpha, step = step,
                            smooth = FALSE, verbose = FALSE)
  return(fit$clusters)
}

cluster_fpcac <- function(Z, K, trim = 0.05, nbasis = 15, nharmonics = 3) {
  # FPCAC (Adelfio et al. 2011): B-spline FPCA + trimmed k-means, as
  # implemented in clustEff::fpcac(). Trimmed curves receive label 0.
  if (!requireNamespace("clustEff", quietly = TRUE)) {
    stop("Package 'clustEff' is required.")
  }
  fit <- clustEff::fpcac(t(Z), K = K, nbasis = nbasis, nharmonics = nharmonics,
                         alpha = trim)
  return(fit$clusters)
}

available_algorithms <- function() {
  # Typing algorithms whose packages are installed on this machine.
  out <- c("fpca_kmeans", "fpca_ward")
  if (requireNamespace("tclust", quietly = TRUE)) out <- c(out, "fpca_tkmeans")
  if (requireNamespace("mclust", quietly = TRUE)) out <- c(out, "fpca_mclust")
  if (requireNamespace("cluster", quietly = TRUE)) out <- c(out, "fpca_pam")
  if (requireNamespace("clustEff", quietly = TRUE)) {
    out <- c(out, "fpcac", "clusteff")
  }
  return(out)
}

choose_K_silhouette <- function(scores, labels_by_K) {
  # labels_by_K: list of label vectors, one per candidate K.
  # Returns the index of the candidate with the largest mean silhouette.
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Package 'cluster' is required for the silhouette.")
  }
  D <- dist(scores)
  mean_sil <- numeric(length(labels_by_K))
  for (j in seq_along(labels_by_K)) {
    lab <- labels_by_K[[j]]
    ok <- lab > 0
    sil <- cluster::silhouette(lab[ok], dist(scores[ok, , drop = FALSE]))
    mean_sil[j] <- mean(sil[, "sil_width"])
  }
  return(list(best = which.max(mean_sil), mean_silhouette = mean_sil))
}

type_curves <- function(Y, K, representation = "spectrum",
                        algorithm = "fpca_kmeans", var_explained = 0.90,
                        trim = 0.05) {
  # Full typing step: representation -> (FPCA) -> clustering -> canonical
  # relabelling. algorithm in {"fpca_kmeans", "fpca_tkmeans", "fpca_mclust",
  # "fpca_ward", "fpca_pam", "ward", "pam", "fpcac", "clusteff"}; "ward" and "pam"
  # work on the full representation, the "fpca_" versions on FPCA scores.
  Z <- represent_curves(Y, method = representation)
  if (algorithm == "fpca_kmeans") {
    S <- fpca_scores(Z, var_explained)
    labels <- cluster_kmeans(S, K)
  } else if (algorithm == "fpca_tkmeans") {
    S <- fpca_scores(Z, var_explained)
    labels <- cluster_trimmed_kmeans(S, K, trim)
  } else if (algorithm == "fpca_mclust") {
    S <- fpca_scores(Z, var_explained)
    labels <- cluster_mclust(S, K)
  } else if (algorithm == "ward") {
    labels <- cluster_ward(Z, K)
  } else if (algorithm == "fpca_ward") {
    labels <- cluster_ward(fpca_scores(Z, var_explained), K)
  } else if (algorithm == "pam") {
    labels <- cluster_pam(Z, K)
  } else if (algorithm == "fpca_pam") {
    labels <- cluster_pam(fpca_scores(Z, var_explained), K)
  } else if (algorithm == "clusteff") {
    labels <- cluster_clusteff(Z, K)
  } else if (algorithm == "fpcac") {
    labels <- cluster_fpcac(Z, K, trim)
  } else {
    stop("Unknown algorithm: ", algorithm)
  }
  labels <- relabel_canonical(labels, Z)
  return(labels)
}

adjusted_rand_index <- function(a, b) {
  # Hubert & Arabie (1985) adjusted Rand index between two partitions.
  tab <- table(a, b)
  comb2 <- function(x) x * (x - 1) / 2
  sum_ij <- sum(comb2(tab))
  sum_a <- sum(comb2(rowSums(tab)))
  sum_b <- sum(comb2(colSums(tab)))
  expected <- sum_a * sum_b / comb2(length(a))
  maximum <- (sum_a + sum_b) / 2
  return((sum_ij - expected) / (maximum - expected))
}


## =============================================================================
## 5. SPATIAL INGREDIENTS
## =============================================================================

cv_bandwidth_types <- function(X, labels, grid = NULL) {
  # Likelihood cross-validation bandwidth for kernel estimates of the type
  # probabilities (as the default of spatstat::bw.relrisk): maximises
  #   sum_i log p_hat_{c_i}(x_i),  p_hat leave-one-out Nadaraya-Watson.
  if (is.null(grid)) grid <- default_bandwidth(X) * exp(seq(log(0.25), log(4), length.out = 15))
  D2 <- as.matrix(dist(X))^2
  same <- outer(labels, labels, "==")
  best <- grid[1]
  best_val <- -Inf
  for (h in grid) {
    W <- exp(-D2 / (2 * h^2))
    diag(W) <- 0
    p_own <- rowSums(W * same) / rowSums(W)
    val <- sum(log(pmax(p_own, 1e-10)))
    if (val > best_val) {
      best_val <- val
      best <- h
    }
  }
  return(best)
}

default_bandwidth <- function(X) {
  # Simple normal-reference bandwidth for a bivariate Gaussian kernel.
  n <- nrow(X)
  s <- mean(c(sd(X[, 1]), sd(X[, 2])))
  return(s * n^(-1 / 6))
}

kde_intensity <- function(X, bw = NULL) {
  # Leave-one-out Gaussian kernel estimate of the intensity at the data
  # points, without edge correction (sufficient for weights; see paper).
  if (is.null(bw)) bw <- default_bandwidth(X)
  D2 <- as.matrix(dist(X))^2
  K <- exp(-D2 / (2 * bw^2)) / (2 * pi * bw^2)
  diag(K) <- 0
  return(rowSums(K))
}

kde_at <- function(points, X, bw = NULL) {
  # Gaussian kernel intensity estimate at arbitrary locations 'points'
  # (m x 2 matrix) from data X. Used only by the centroid statistic.
  if (is.null(bw)) bw <- default_bandwidth(X)
  points <- matrix(points, ncol = 2)
  out <- numeric(nrow(points))
  for (m in 1:nrow(points)) {
    d2 <- (X[, 1] - points[m, 1])^2 + (X[, 2] - points[m, 2])^2
    out[m] <- sum(exp(-d2 / (2 * bw^2))) / (2 * pi * bw^2)
  }
  return(out)
}

default_scales <- function(X, kernel = "ball") {
  # Default spatial scales, relative to the side of the bounding square:
  #   ball kernel     : radii 0.05, 0.10, 0.15, 0.20 times the side
  #   gaussian kernel : bandwidths 1/32, 1/16, 1/8 times the side
  side <- sqrt(diff(range(X[, 1])) * diff(range(X[, 2])))
  if (kernel == "ball") return(side * c(0.05, 0.10, 0.15, 0.20))
  return(side * c(1, 2, 4) / 32)
}

build_pair_list <- function(X, scales, weights = "none", bw = NULL,
                            kernel = "ball") {
  # Pre-computes everything that does NOT depend on the labels:
  #   - the list of unordered pairs (i < j) that can receive a weight,
  #   - the pair weights a_ij(s) for every scale s (matrix |pairs| x S):
  #       kernel = "ball"     : a_ij(r) = w_ij * 1{d_ij <= r}       (scale = radius)
  #       kernel = "gaussian" : a_ij(h) = w_ij * exp(-d_ij^2 / (2 h^2))
  #                                                             (scale = bandwidth;
  #                             pairs farther than 4 h_max get weight ~ 0 and
  #                             are dropped)
  #   - the quantities S0, S1, S2 of the closed-form null moments, per scale.
  # weights: "none" (w_ij = 1) or "inverse_intensity"
  #          (w_ij = 1 / (lambda_i * lambda_j), rescaled to mean 1).
  # The element 'radii' (= scales) and, for the ball kernel, 'bin' are kept
  # for backward compatibility.
  n <- nrow(X)
  scales <- sort(scales)
  R <- length(scales)
  D <- as.matrix(dist(X))
  if (kernel == "ball") {
    reach <- max(scales)
  } else if (kernel == "gaussian") {
    reach <- 4 * max(scales)
  } else {
    stop("Unknown kernel: ", kernel)
  }
  idx <- which(upper.tri(D) & D <= reach, arr.ind = TRUE)
  i <- idx[, 1]
  j <- idx[, 2]
  d <- D[idx]
  if (weights == "none") {
    w <- rep(1, length(d))
    lambda <- rep(NA_real_, n)
  } else if (weights == "inverse_intensity") {
    lambda <- kde_intensity(X, bw)
    w <- 1 / (lambda[i] * lambda[j])
    w <- w / mean(w)
  } else {
    stop("Unknown weights: ", weights)
  }
  a <- matrix(0, nrow = length(d), ncol = R)
  for (m in 1:R) {
    if (kernel == "ball") {
      a[, m] <- w * (d <= scales[m])
    } else {
      a[, m] <- w * exp(-d^2 / (2 * scales[m]^2))
    }
  }
  bin <- NULL
  if (kernel == "ball") {
    bin <- findInterval(d, c(0, scales), left.open = TRUE, rightmost.closed = TRUE)
    bin[d == 0] <- 1
  }
  # Null-moment ingredients for each scale (ordered-pair convention):
  #   S0 = sum_{i != j} a_ij, S1 = sum_{i != j} a_ij^2,
  #   S2 = sum_i (sum_{j != i} a_ij)^2.
  S0 <- numeric(R)
  S1 <- numeric(R)
  S2 <- numeric(R)
  for (m in 1:R) {
    S0[m] <- 2 * sum(a[, m])
    S1[m] <- 2 * sum(a[, m]^2)
    row_sums <- tabulate_weighted(i, a[, m], n) + tabulate_weighted(j, a[, m], n)
    S2[m] <- sum(row_sums^2)
  }
  return(list(n = n, i = i, j = j, d = d, w = w, a = a, bin = bin,
              radii = scales, scales = scales, kernel = kernel,
              S0 = S0, S1 = S1, S2 = S2, lambda = lambda, X = X))
}

tabulate_weighted <- function(index, weight, n) {
  # Sum of 'weight' by integer 'index' in 1..n (weighted tabulate).
  out <- numeric(n)
  if (length(index) == 0) return(out)
  s <- rowsum(weight, index)
  out[as.integer(rownames(s))] <- s[, 1]
  return(out)
}


## =============================================================================
## 6. TEST STATISTICS
## =============================================================================

pair_statistic <- function(labels, pairs, K) {
  # U_k(s) for k = 1..K and all scales s, ordered-pair convention:
  #   U_k(s) = sum_{i != j} 1{c_i = c_j = k} a_ij(s) = 2 * sum_{i < j} ... .
  # Returns a K x S matrix.
  R <- ncol(pairs$a)
  li <- labels[pairs$i]
  lj <- labels[pairs$j]
  same <- li == lj & li > 0
  U <- matrix(0, nrow = K, ncol = R)
  for (k in 1:K) {
    sel <- same & li == k
    if (any(sel)) {
      U[k, ] <- 2 * colSums(pairs$a[sel, , drop = FALSE])
    }
  }
  return(U)
}

pair_null_moments <- function(n, n_k, pairs) {
  # Exact mean and variance of U_k(r) under uniform random permutation of
  # the labels (n_k points of type k among n), Lemma S1 of the Supplementary Material:
  #   E U   = p2 * S0
  #   E U^2 = 2 p2 S1 + 4 p3 (S2 - S1) + p4 (S0^2 + 2 S1 - 4 S2)
  # with p_m = n_k (n_k - 1) ... (n_k - m + 1) / [n (n - 1) ... (n - m + 1)].
  p2 <- n_k * (n_k - 1) / (n * (n - 1))
  p3 <- p2 * (n_k - 2) / (n - 2)
  p4 <- p3 * (n_k - 3) / (n - 3)
  p3 <- max(p3, 0)
  p4 <- max(p4, 0)
  mean_U <- p2 * pairs$S0
  second <- 2 * p2 * pairs$S1 + 4 * p3 * (pairs$S2 - pairs$S1) +
    p4 * (pairs$S0^2 + 2 * pairs$S1 - 4 * pairs$S2)
  var_U <- pmax(second - mean_U^2, 0)
  return(list(mean = mean_U, var = var_U))
}

all_null_moments <- function(labels, pairs, K) {
  # Null moments for every type k. They depend on the labels only through
  # n_k, which is invariant under permutation: compute them once.
  # Note: elements stay NULL for types with fewer than 2 points.
  moments <- vector("list", K)
  for (k in 1:K) {
    n_k <- sum(labels == k)
    if (n_k >= 2) {
      moments[[k]] <- pair_null_moments(pairs$n, n_k, pairs)
    }
  }
  return(moments)
}

standardise_pair <- function(U, moments, K) {
  # Z_k(r) = (U_k(r) - E_0 U_k(r)) / sd_0 U_k(r).
  # Standardisation only puts radii and types on a common scale; validity
  # comes from the permutation calibration, also under restricted
  # (stratified) permutations whose moments differ from these.
  Zs <- matrix(NA_real_, nrow = K, ncol = ncol(U))
  for (k in 1:K) {
    if (!is.null(moments[[k]])) {
      mom <- moments[[k]]
      ok <- mom$var > 0
      Zs[k, ok] <- (U[k, ok] - mom$mean[ok]) / sqrt(mom$var[ok])
    }
  }
  return(Zs)
}

coordinatewise_median <- function(P) {
  return(c(median(P[, 1]), median(P[, 2])))
}

centroid_statistic <- function(labels, X, K, radii, lambda_points = NULL,
                               bw = NULL) {
  # Centroid-based dispersion statistic (Section 3.1 of the paper):
  #   S(k) = sum_j sum_r 1{d(M_k, x_j) > r} 1{x_j in C_k} /
  #          (|C_k| lambda(M_k) lambda(x_j)),
  # M_k = coordinate-wise median of cluster k, lambda = Gaussian KDE.
  # Small values = spatially compact cluster (left-tailed test).
  if (is.null(bw)) bw <- default_bandwidth(X)
  if (is.null(lambda_points)) lambda_points <- kde_at(X, X, bw)
  S <- rep(NA_real_, K)
  for (k in 1:K) {
    members <- which(labels == k)
    if (length(members) >= 2) {
      M <- coordinatewise_median(X[members, , drop = FALSE])
      lambda_M <- kde_at(M, X, bw)
      d <- sqrt((X[members, 1] - M[1])^2 + (X[members, 2] - M[2])^2)
      total <- 0
      for (r in radii) {
        total <- total + sum((d > r) / lambda_points[members])
      }
      S[k] <- total / (length(members) * lambda_M)
    }
  }
  return(S)
}


## =============================================================================
## 7. PERMUTATION MACHINERY
## =============================================================================

draw_permutation <- function(n, strata = NULL) {
  # Uniform permutation of 1..n; if 'strata' is given (vector of length n),
  # indices are permuted only within each stratum (restricted
  # randomisation for the covariate-adjusted null, Proposition 1).
  if (is.null(strata)) {
    return(sample.int(n))
  }
  perm <- 1:n
  for (s in unique(strata)) {
    members <- which(strata == s)
    if (length(members) > 1) {
      perm[members] <- members[sample.int(length(members))]
    }
  }
  return(perm)
}

run_replicates <- function(one_replicate, n_rep, setting_id, n_cores = 1,
                           seed_base = 0) {
  # Runs one_replicate() n_rep times and returns a list of data frames.
  # Replicate r of setting s uses the seed seed_base + 100000 * s + r, set
  # INSIDE the replicate: results are reproducible whatever the number of
  # cores, and replicates of different settings are independent.
  #
  # Robustness: an error inside a replicate is caught and reported with the
  # replicate number; replicates that fail in parallel (error or crashed
  # worker) are rerun serially, so that a real error stops the script with
  # an informative message instead of producing empty results.
  worker <- function(r) {
    out <- tryCatch({
      set.seed(seed_base + 100000 * setting_id + r)
      row <- one_replicate()
      row$rep <- r
      row
    }, error = function(e) {
      paste0("replicate ", r, ": ", conditionMessage(e))
    })
    return(out)
  }
  if (n_cores > 1) {
    rows <- parallel::mclapply(1:n_rep, worker, mc.cores = n_cores,
                               mc.preschedule = TRUE)
  } else {
    rows <- list()
    for (r in 1:n_rep) rows[[r]] <- worker(r)
  }
  failed <- which(!vapply(rows, is.data.frame, logical(1)))
  if (length(failed) > 0) {
    first_msg <- rows[[failed[1]]]
    if (!is.character(first_msg)) first_msg <- "worker process did not return a result"
    message("setting ", setting_id, ": ", length(failed),
            " replicate(s) failed (first: ", first_msg, "); rerunning serially")
    for (r in failed) {
      rows[[r]] <- worker(r)
      if (!is.data.frame(rows[[r]])) {
        stop("setting ", setting_id, ", ", rows[[r]], call. = FALSE)
      }
    }
  }
  return(rows)
}

make_strata <- function(z, n_bins = 5) {
  # Quantile bins of a numeric covariate (or interaction of several).
  if (is.null(dim(z))) z <- matrix(z, ncol = 1)
  codes <- rep("", nrow(z))
  for (col in 1:ncol(z)) {
    br <- unique(quantile(z[, col], probs = seq(0, 1, length.out = n_bins + 1)))
    b <- cut(z[, col], breaks = br, include.lowest = TRUE, labels = FALSE)
    codes <- paste(codes, b, sep = "_")
  }
  return(as.integer(factor(codes)))
}


## =============================================================================
## 8. THE TEST
## =============================================================================

holm_adjust <- function(p) {
  # Holm (1979) step-down adjustment; NA entries are ignored.
  out <- rep(NA_real_, length(p))
  ok <- !is.na(p)
  out[ok] <- p.adjust(p[ok], method = "holm")
  return(out)
}

format_seconds <- function(sec) {
  # 3725 -> "1h 02m", 125 -> "2m 05s", 42 -> "42s"
  sec <- round(sec)
  if (sec >= 3600) return(sprintf("%dh %02dm", sec %/% 3600, (sec %% 3600) %/% 60))
  if (sec >= 60) return(sprintf("%dm %02ds", sec %/% 60, sec %% 60))
  return(sprintf("%ds", sec))
}

progress_report <- function(b, B, t0, label = "") {
  # Prints, on one line updated in place, the percentage of permutations done,
  # the elapsed time and an estimate of the remaining time (about 20 updates).
  step <- max(1, floor(B / 20))
  if (b %% step == 0 || b == B) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    remaining <- elapsed / b * (B - b)
    cat(sprintf("\r   %s %3.0f%%  (%d/%d)  elapsed %s  remaining ~%s     ",
                label, 100 * b / B, b, B, format_seconds(elapsed),
                format_seconds(remaining)))
    if (b == B) cat("\n")
    flush.console()
  }
}

fcst_test <- function(X, labels, scales = NULL, B = 999, kernel = "ball",
                      weights = "none", strata = NULL, centroid = FALSE,
                      centroid_radii = NULL, bw = NULL, keep_null = FALSE,
                      progress = FALSE, progress_label = "FCST") {
  # Functional Clustering Spatial Test (Section 3 of the paper).
  #
  # Arguments
  #   X        n x 2 coordinates
  #   labels   vector of n labels (0 = unassigned) OR an n x A matrix with
  #            one column per typing (different algorithms / K / features).
  #            All columns are permuted with the SAME permutation, which makes
  #            the max-type calibration exact over the whole family
  #            (Proposition 1).
  #   scales   radii (kernel = "ball", default) or bandwidths
  #            (kernel = "gaussian"); default: default_scales(X, kernel)
  #   B        number of random permutations
  #   kernel   "ball" (default, indicator 1{d <= r}) or "gaussian"
  #   weights  "none" or "inverse_intensity"
  #   strata   NULL or vector of strata for restricted permutations
  #   centroid if TRUE also computes the centroid-based S(k) test on the
  #            first typing, with radii centroid_radii (default
  #            c(0.05, 0.10, 0.15, 0.20) times the side of the window)
  #
  # Statistics, for typing a, type k and scale s:
  #   Z_{a,k}(s) = standardised same-type weight (exact permutation moments)
  #   T_{a,k}    = max_s Z_{a,k}(s)                    type-specific
  #   G_a        = max_s sum_k Z_{a,k}(s)               global, typing a
  #   G          = max_a G_a                            global, all typings
  # p-values, all of the form (1 + #{null >= obs}) / (B + 1):
  #   p_global : test of random labelling based on G (Propositions 1-2)
  #   p_typing : test based on G_a, for each typing separately
  #   table$p_raw  : type-specific test based on T_{a,k}; valid under the
  #                  type-specific null H_0k (Proposition 2)
  #   table$p_holm : Holm adjustment of p_raw across the types of typing a
  #                  (strong control of the FWER over types)
  #   table$p_maxT : single-step max-T adjustment over all (a, k)
  labels <- as.matrix(labels)
  n <- nrow(X)
  A <- ncol(labels)
  if (is.null(scales)) scales <- default_scales(X, kernel)
  pairs <- build_pair_list(X, scales, weights = weights, bw = bw,
                           kernel = kernel)
  K_by <- integer(A)
  moments_by <- list()
  for (a in 1:A) {
    K_by[a] <- max(labels[, a])
    moments_by[[a]] <- all_null_moments(labels[, a], pairs, K_by[a])
  }

  # ---- statistics for a given label matrix -------------------------------------
  # returns c(T_{1,1..K1}, ..., T_{A,1..KA}, G_1, ..., G_A)
  compute_stats <- function(lab_matrix) {
    T_all <- numeric(0)
    G_all <- numeric(A)
    for (a in 1:A) {
      U <- pair_statistic(lab_matrix[, a], pairs, K_by[a])
      Zs <- standardise_pair(U, moments_by[[a]], K_by[a])
      Tk <- rep(NA_real_, K_by[a])
      for (k in 1:K_by[a]) {
        if (any(!is.na(Zs[k, ]))) Tk[k] <- max(Zs[k, ], na.rm = TRUE)
      }
      T_all <- c(T_all, Tk)
      G_all[a] <- max(colSums(Zs, na.rm = TRUE))
    }
    return(c(T_all, G_all))
  }

  obs <- compute_stats(labels)
  n_type_stats <- sum(K_by)
  T_obs <- obs[1:n_type_stats]
  G_obs <- obs[n_type_stats + (1:A)]
  T_null <- matrix(NA_real_, nrow = B, ncol = n_type_stats)
  G_null <- matrix(NA_real_, nrow = B, ncol = A)

  if (centroid) {
    if (is.null(centroid_radii)) {
      side <- sqrt(diff(range(X[, 1])) * diff(range(X[, 2])))
      centroid_radii <- side * c(0.05, 0.10, 0.15, 0.20)
    }
    lambda_points <- kde_at(X, X, bw)
    S_obs <- centroid_statistic(labels[, 1], X, K_by[1], centroid_radii,
                                lambda_points, bw)
    S_null <- matrix(NA_real_, nrow = B, ncol = K_by[1])
  }

  t0 <- Sys.time()
  for (b in 1:B) {
    if (progress) progress_report(b, B, t0, progress_label)
    perm <- draw_permutation(n, strata)
    lab_b <- labels[perm, , drop = FALSE]
    stats_b <- compute_stats(lab_b)
    T_null[b, ] <- stats_b[1:n_type_stats]
    G_null[b, ] <- stats_b[n_type_stats + (1:A)]
    if (centroid) {
      S_null[b, ] <- centroid_statistic(lab_b[, 1], X, K_by[1], centroid_radii,
                                        lambda_points, bw)
    }
  }

  # ---- global p-values -------------------------------------------------------------
  p_typing <- numeric(A)
  for (a in 1:A) {
    p_typing[a] <- (1 + sum(G_null[, a] >= G_obs[a])) / (B + 1)
  }
  G_max_null <- apply(G_null, 1, max)
  p_global <- (1 + sum(G_max_null >= max(G_obs))) / (B + 1)

  # ---- type-specific p-values ------------------------------------------------------
  p_raw <- rep(NA_real_, n_type_stats)
  p_maxT <- rep(NA_real_, n_type_stats)
  max_null <- apply(T_null, 1, max, na.rm = TRUE)
  for (s in 1:n_type_stats) {
    if (!is.na(T_obs[s])) {
      p_raw[s] <- (1 + sum(T_null[, s] >= T_obs[s], na.rm = TRUE)) / (B + 1)
      p_maxT[s] <- (1 + sum(max_null >= T_obs[s])) / (B + 1)
    }
  }
  typing <- integer(0)
  type <- integer(0)
  for (a in 1:A) {
    typing <- c(typing, rep(a, K_by[a]))
    type <- c(type, 1:K_by[a])
  }
  p_holm <- rep(NA_real_, n_type_stats)
  for (a in 1:A) {
    idx <- which(typing == a)
    p_holm[idx] <- holm_adjust(p_raw[idx])
  }
  table_out <- data.frame(typing = typing, type = type,
                          n_k = NA_integer_, T_obs = T_obs,
                          p_raw = p_raw, p_holm = p_holm, p_maxT = p_maxT)
  for (s in 1:n_type_stats) {
    table_out$n_k[s] <- sum(labels[, typing[s]] == type[s])
  }

  result <- list(table = table_out, p_global = p_global,
                 p_typing = p_typing, G_obs = G_obs, B = B,
                 scales = scales, kernel = kernel, weights = weights,
                 stratified = !is.null(strata))
  if (centroid) {
    p_centroid <- rep(NA_real_, K_by[1])
    for (k in 1:K_by[1]) {
      if (!is.na(S_obs[k])) {
        # left-tailed: small S(k) = compact cluster
        p_centroid[k] <- (1 + sum(S_null[, k] <= S_obs[k], na.rm = TRUE)) / (B + 1)
      }
    }
    result$centroid <- data.frame(type = 1:K_by[1], S_obs = S_obs,
                                  p_raw = p_centroid,
                                  p_bonf = pmin(1, K_by[1] * p_centroid))
  }
  if (keep_null) {
    result$T_null <- T_null
    result$G_null <- G_null
  }
  return(result)
}


## =============================================================================
## 9. MISCLASSIFICATION HELPERS (Proposition 3)
## =============================================================================

add_label_noise <- function(labels, K, eps) {
  # Symmetric non-differential misclassification: each label is replaced,
  # with probability eps, by one of the other K - 1 labels chosen uniformly.
  n <- length(labels)
  noisy <- labels
  flip <- runif(n) < eps
  for (i in which(flip)) {
    others <- setdiff(1:K, labels[i])
    noisy[i] <- others[sample.int(K - 1, 1)]
  }
  return(noisy)
}

attenuation_gamma <- function(eps, K) {
  # p_hat_k(x) = eps / (K - 1) + gamma * p_k(x),  gamma = 1 - eps K / (K - 1)
  return(1 - eps * K / (K - 1))
}

delta_signal <- function(p_k, pairs, scale_index) {
  # Population signal functional (Lemma 1 in the paper), for any kernel:
  #   Delta_k(s) = sum_{i != j} a_ij(s) (q_i - q_bar)(q_j - q_bar),
  #   q_bar = sum_i a_i. q_i / sum_i a_i.,  q_i = p_k(x_i).
  i <- pairs$i
  j <- pairs$j
  a <- pairs$a[, scale_index]
  n <- pairs$n
  a_dot <- tabulate_weighted(i, a, n) + tabulate_weighted(j, a, n)
  q_bar <- sum(a_dot * p_k) / sum(a_dot)
  dev <- p_k - q_bar
  return(2 * sum(a * dev[i] * dev[j]))
}


## =============================================================================
## 10. LOCAL DIAGNOSTICS (descriptive localisation of hotspots)
## =============================================================================

local_same_type <- function(labels, pairs, scale_index, B = 999) {
  # For each point i: L_i = sum_{j != i} a_ij(s) 1{c_j = c_i}.
  # Conditional permutation p-value: c_i is held fixed and the other labels
  # are permuted (as in LISA-type statistics). Returns L and p-values; the
  # p-values are point-wise and must be corrected (e.g. Benjamini-Hochberg)
  # before mapping.
  n <- pairs$n
  in_r <- pairs$a[, scale_index] > 0
  i <- pairs$i[in_r]
  j <- pairs$j[in_r]
  a <- pairs$a[in_r, scale_index]

  local_counts <- function(lab) {
    same <- lab[i] == lab[j] & lab[i] > 0
    tabulate_weighted(i[same], a[same], n) + tabulate_weighted(j[same], a[same], n)
  }

  L_obs <- local_counts(labels)
  exceed <- rep(0, n)
  for (b in 1:B) {
    # Approximation of the conditional scheme (documented in the paper):
    # one global permutation per replicate; when point pt is evaluated its
    # own label is kept fixed and the neighbours take the permuted labels.
    perm <- sample.int(n)
    lab_b <- labels[perm]
    L_b <- numeric(n)
    for (pt in 1:n) {
      nb_left <- j[i == pt]
      nb_right <- i[j == pt]
      nb <- c(nb_left, nb_right)
      wt <- c(a[i == pt], a[j == pt])
      L_b[pt] <- sum(wt[lab_b[nb] == labels[pt]])
    }
    exceed <- exceed + (L_b >= L_obs)
  }
  p_local <- (1 + exceed) / (B + 1)
  return(data.frame(L = L_obs, p_local = p_local,
                    p_bh = p.adjust(p_local, method = "BH")))
}
