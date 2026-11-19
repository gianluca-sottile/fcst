################################################################################
#  make_design_figures.R
#
#  Illustrative figures of the simulation design (Section 6.1 of the paper):
#    fig_design_space.pdf  : one realisation of each spatial setting
#                            (top: point processes under random labelling;
#                             bottom: the three alternatives, beta = 4)
#    fig_design_curves.pdf : example curves of the three types at the three
#                            signal-to-noise ratios used in the simulations,
#                            with the amplitude- and phase-invariant
#                            representation (Fourier modulus) underneath
#    fig_confounding.pdf   : one realisation of the confounding scenario
#    fig_design_robust.pdf : hotspot profiles and settings of sim_08
#                            (Supplementary Material)
#
#  All figures are drawn at their final printed size with the common style of
#  plot_style.R (9-point font, Okabe-Ito colours).
#
#  Base R only. Run from the R/ folder:  source("make_design_figures.R")
################################################################################

source("fcst_methods.R")
source("plot_style.R")
dir_fig <- file.path("..", "output", "figures")
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

set.seed(20261030)
K <- 3

type_cols   <- col_types[1:K]
type_pch    <- pch_types[1:K]
cex_type    <- c(0.75, 0.65, 0.6)        # type 1 slightly larger
col_station <- "black"
type_labels <- c("type 1", "type 2", "type 3")

## ---- helper: one point pattern coloured by type ----------------------------------
plot_pattern <- function(X, types, letter, main, surface = NULL) {
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), asp = 1, xlab = "", ylab = "",
       xaxt = "n", yaxt = "n", bty = "n")
  panel_title(letter, main)
  if (!is.null(surface)) {
    # contours of the type-1 excess surface h_1(x)
    g <- seq(0, 1, length.out = 151)
    grid_pts <- as.matrix(expand.grid(g, g))
    h <- matrix(hotspot_surface(grid_pts, surface), nrow = 151)
    contour(g, g, h, levels = c(0.25, 0.75), add = TRUE, drawlabels = FALSE,
            col = col_reference, lty = 2, lwd = 1)
  }
  for (k in c(3, 2, 1)) {            # type 1 drawn last, on top
    sel <- types == k
    points(X[sel, 1], X[sel, 2], col = type_cols[k], pch = type_pch[k],
           cex = cex_type[k])
  }
  rect(0, 0, 1, 1, lwd = lwd_box)
}

## =============================================================================
## Figure: spatial settings
## =============================================================================
n <- 300
open_pdf(file.path(dir_fig, "fig_design_space.pdf"), height = 4.75)
set_panels(c(2, 3), mar = c(0.3, 0.3, 1.5, 0.3), oma = c(1.6, 0, 0, 0))

process_titles <- c(hpp = "Homogeneous Poisson", ipp = "Inhomogeneous Poisson",
                    thomas = "Thomas process")
panel <- 0
for (proc in c("hpp", "ipp", "thomas")) {
  panel <- panel + 1
  X <- sim_points(n, proc)
  types <- sim_types(type_probabilities(X, K = K, beta = 0))
  plot_pattern(X, types, letters[panel], process_titles[[proc]])
}
pattern_titles <- c(one_bump = "One hotspot", two_bumps = "Two hotspots",
                    band = "Elongated band")
for (pat in c("one_bump", "two_bumps", "band")) {
  panel <- panel + 1
  X <- sim_points(n, "hpp")
  types <- sim_types(type_probabilities(X, K = K, beta = 4, pattern = pat))
  plot_pattern(X, types, letters[panel], pattern_titles[[pat]], surface = pat)
}
legend_bottom(legend = c(type_labels, expression("contours " * h[1] == 0.25 * ", " * 0.75)),
              col = c(type_cols, col_reference), pch = c(type_pch, NA),
              pt.cex = c(1, 0.9, 0.9, NA), lty = c(NA, NA, NA, 2),
              lwd = c(NA, NA, NA, 1), ncol = 4)
dev.off()
cat("written:", file.path(dir_fig, "fig_design_space.pdf"), "\n")

## =============================================================================
## Figure: curves and their invariant representation
## =============================================================================
## Top row: two noisy curves per type (thin, semi-transparent), offset
## vertically by type, with the noiseless signal of unit amplitude and zero
## phase (thick). Bottom row: Fourier modulus of 30 curves per type (thin)
## and the type mean (thick).
snr_levels <- c(1, 0.05, 0.03)
shape <- type_shape_parameters(K)

open_pdf(file.path(dir_fig, "fig_design_curves.pdf"), height = 4.6)
set_panels(c(2, 3), mar = c(2.8, 3.4, 1.5, 0.4), oma = c(1.6, 0, 0, 0))
panel <- 0
for (s in seq_along(snr_levels)) {
  panel <- panel + 1
  types <- rep(1:K, each = 2)
  Y <- sim_curves(types, snr = snr_levels[s])
  tt <- attr(Y, "grid")
  step <- 4 * sd(as.vector(Y))
  offset <- c(0, -1, -2) * step
  plot(NA, xlim = c(0, 1), ylim = c(-2.6 * step, 0.6 * step),
       xlab = "t", ylab = "", yaxt = "n", bty = "n")
  panel_title(letters[panel], paste0("Curves, SNR ", snr_levels[s]))
  for (i in seq_along(types)) {
    lines(tt, Y[i, ] + offset[types[i]],
          col = transparent(type_cols[types[i]], 0.45), lwd = lwd_thin)
  }
  for (k in 1:K) {
    envelope <- (tt / shape$peak[k]) * exp(1 - tt / shape$peak[k])
    signal <- envelope * sin(2 * pi * shape$freq[k] * tt)
    lines(tt, offset[k] + signal, col = type_cols[k], lwd = lwd_main)
  }
  axis(2, at = offset, labels = paste("type", 1:K), tick = FALSE)
  box(lwd = lwd_box)
}
for (s in seq_along(snr_levels)) {
  panel <- panel + 1
  types <- rep(1:K, each = 30)
  Y <- sim_curves(types, snr = snr_levels[s])
  Z <- represent_curves(Y, method = "spectrum")
  freq <- 1:ncol(Z)
  plot(NA, xlim = c(1, 25), ylim = c(0, max(Z[, 1:25])), bty = "n",
       xlab = "frequency (cycles on [0, 1])", ylab = "normalised modulus")
  panel_title(letters[panel], paste0("Spectra, SNR ", snr_levels[s]))
  for (i in seq_along(types)) {
    lines(freq, Z[i, ], col = transparent(type_cols[types[i]], 0.15),
          lwd = lwd_thin)
  }
  for (k in 1:K) {
    lines(freq, colMeans(Z[types == k, , drop = FALSE]), col = type_cols[k],
          lwd = lwd_main)
  }
  box(lwd = lwd_box)
}
legend_bottom(legend = c("type 1", "type 2", "type 3"),
              col = type_cols, lwd = lwd_main, lty = 1, ncol = 3)
dev.off()
cat("written:", file.path(dir_fig, "fig_design_curves.pdf"), "\n")

## =============================================================================
## Figure: confounding scenario
## =============================================================================
## (a) true types, independent of location; (b) clusters obtained from the
## raw curves: the minority clusters collect a few high-amplitude events
## close to the station; (c) curve norm versus distance from the station,
## coloured as in (b).
n <- 300
X <- sim_points(n, "hpp")
z <- sqrt((X[, 1] - 0.5)^2 + (X[, 2] - 0.5)^2)
types <- sim_types(type_probabilities(X, K = K, beta = 0))
Y <- sim_curves(types, snr = 1, amp_multiplier = exp(-4 * z))
lab_raw  <- type_curves(Y, K, representation = "raw", algorithm = "fpca_kmeans")
lab_spec <- type_curves(Y, K, representation = "spectrum",
                        algorithm = "fpca_kmeans")
curve_norm <- sqrt(rowSums(Y^2))
fmt_ari <- function(x) {
  x <- round(x, 2)
  if (x == 0) x <- 0          # avoid printing "-0.00"
  formatC(x, format = "f", digits = 2)
}
ari_raw  <- fmt_ari(adjusted_rand_index(types, lab_raw))
ari_spec <- fmt_ari(adjusted_rand_index(types, lab_spec))
cat("confounding example: ARI raw =", ari_raw, " ARI spectrum =", ari_spec, "\n")

majority <- as.integer(names(which.max(table(lab_raw))))
is_minor <- lab_raw != majority
col_major <- "grey70"
col_minor <- type_cols[1]

open_pdf(file.path(dir_fig, "fig_confounding.pdf"), height = 2.7)
set_panels(c(1, 3), mar = c(0.3, 0.3, 1.5, 0.3), oma = c(2.6, 0, 0, 0))
plot_pattern(X, types, "a", "True types")
points(0.5, 0.5, pch = 4, cex = 1.6, lwd = 2.2, col = col_station)

plot(NA, xlim = c(0, 1), ylim = c(0, 1), asp = 1, xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n")
panel_title("b", paste0("Raw-curve clusters, ARI ", ari_raw))
points(X[!is_minor, 1], X[!is_minor, 2], col = col_major, pch = 16, cex = 0.6)
points(X[is_minor, 1], X[is_minor, 2], col = col_minor, pch = 16, cex = 1.1)
points(0.5, 0.5, pch = 4, cex = 1.6, lwd = 2.2, col = col_station)
rect(0, 0, 1, 1, lwd = lwd_box)

par(mar = c(2.9, 3.3, 1.5, 0.4))
plot(z, curve_norm, log = "y", type = "n", bty = "n",
     xlab = "distance from the station", ylab = "curve norm")
panel_title("c", "Amplitude and distance")
points(z[!is_minor], curve_norm[!is_minor], col = col_major, pch = 16, cex = 0.6)
points(z[is_minor], curve_norm[is_minor], col = col_minor, pch = 16, cex = 1.1)
box(lwd = lwd_box)

legend_bottom(legend = c("type 1", "type 2", "type 3", "station",
                         "majority cluster", "minority clusters"),
              col = c(type_cols, col_station, col_major, col_minor),
              pch = c(type_pch, 4, 16, 16),
              pt.cex = c(1, 0.9, 0.9, 1.3, 0.9, 1.2),
              pt.lwd = c(1, 1, 1, 2, 1, 1), ncol = 3)
dev.off()
cat("written:", file.path(dir_fig, "fig_confounding.pdf"), "\n")

## =============================================================================
## Figure: settings of the robustness study (sim_08, Supplementary Material)
## =============================================================================
## Own seed, so that the figures above are unchanged.
set.seed(20261111)
open_pdf(file.path(dir_fig, "fig_design_robust.pdf"), height = 6.6)
set_panels(c(2, 2), mar = c(0.3, 0.3, 1.5, 0.3), oma = c(1.6, 0, 0, 0))

par(mar = c(3, 3.3, 1.5, 0.6))
s_bump <- 0.10
d_grid <- seq(0, 0.4, length.out = 400)
pts <- cbind(0.3 + d_grid, rep(0.3, length(d_grid)))
h_g <- bump_profile(pts, c(0.3, 0.3), s_bump, "gauss")
h_t <- bump_profile(pts, c(0.3, 0.3), s_bump, "t")
plot(d_grid, h_t, type = "l", lwd = lwd_main, col = type_cols[1], bty = "n",
     xlab = "distance from the centre", ylab = expression(h[1]), ylim = c(0, 3.05))
panel_title("a", "Hotspot profiles")
lines(d_grid, h_g, lwd = lwd_main, col = "grey30", lty = 2)
legend("topright", bty = "n", legend = c("Student t", "Gaussian"),
       col = c(type_cols[1], "grey30"), lty = c(1, 2), lwd = lwd_main)
box(lwd = lwd_box)

par(mar = c(0.3, 0.3, 1.5, 0.3))
X <- sim_points(n, "hpp")
types <- sim_types(type_probabilities(X, K = K, beta = 4, pattern = "one_bump_t"))
plot_pattern(X, types, "b", "Student-t hotspot, Poisson", surface = "one_bump_t")

X <- sim_points(n, "ipp")
types <- sim_types(type_probabilities(X, K = K, beta = 4, pattern = "one_bump"))
plot_pattern(X, types, "c", "Inhomogeneous Poisson", surface = "one_bump")

X <- sim_points_thomas(n, n_parents = 15, sigma = 0.03,
                       fixed_parents = hotspot_centres("two_bumps"))
types <- sim_types(type_probabilities(X, K = K, beta = 4, pattern = "two_bumps"))
plot_pattern(X, types, "d", "Thomas, parents at the hotspots", surface = "two_bumps")

legend_bottom(legend = c(type_labels, expression("contours " * h[1] == 0.25 * ", " * 0.75)),
              col = c(type_cols, col_reference), pch = c(type_pch, NA),
              pt.cex = c(1, 0.9, 0.9, NA), lty = c(NA, NA, NA, 2),
              lwd = c(NA, NA, NA, 1), ncol = 4)
dev.off()
cat("written:", file.path(dir_fig, "fig_design_robust.pdf"), "\n")
