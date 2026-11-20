################################################################################
#  app_05_figures_stead.R
#
#  Figures and descriptive numbers for the results of the STEAD application
#  (Section 7.3 of the paper). Fast: it uses the clusterings saved by
#  app_02_fcst_stead.R (results/app_stead_labels.rds) and does not repeat the
#  tests.
#
#  Output
#    ../output/figures/fig_app_map.pdf     epicentres of each type of the primary
#                                         clustering + zoom on the most
#                                         concentrated type
#    ../output/figures/fig_app_curves.pdf  mean spectra and example waveforms
#    results/app_stead_type_locations.csv where each type is (printed too)
################################################################################

source("fcst_methods.R")
source("plot_style.R")

stead <- readRDS("results/app_stead_data.rds")
lab   <- readRDS("results/app_stead_labels.rds")
types <- read.csv("results/app_stead_types.csv")
X  <- stead$X
Y  <- stead$Y
m  <- stead$meta
lab1 <- lab$labels[, lab$primary]
K1 <- max(lab1)
Z  <- represent_curves(Y, method = "spectrum")
freq <- (1:ncol(Z)) * stead$fs / ncol(Y)

type_cols <- col_types
type_pch  <- pch_types
bg_col    <- col_grey_bg

volcanoes <- data.frame(
  name = c("St. Helens", "Hood", "Adams", "Rainier"),
  lat  = c(46.1912, 45.3736, 46.2024, 46.8528),
  lon  = c(-122.1944, -121.6960, -121.4909, -121.7604))
volcanoes$x <- (volcanoes$lon - stead$lon0) * stead$km_per_deg["lon"]
volcanoes$y <- (volcanoes$lat - stead$lat0) * stead$km_per_deg["lat"]
stations_xy <- unique(stead$station_xy)

## ---- 1. where is each type? --------------------------------------------------------
loc <- data.frame(type = 1:K1, n = as.integer(table(factor(lab1, levels = 1:K1))))
for (v in seq_len(nrow(volcanoes))) {
  dv <- sqrt((X[, 1] - volcanoes$x[v])^2 + (X[, 2] - volcanoes$y[v])^2)
  loc[[paste0("share_within_15km_", gsub("[^A-Za-z]", "", volcanoes$name[v]))]] <-
    as.numeric(tapply(dv <= 15, factor(lab1, levels = 1:K1), mean))
}
loc$median_lon <- as.numeric(tapply(m$source_longitude, factor(lab1, levels = 1:K1), median))
loc$median_lat <- as.numeric(tapply(m$source_latitude,  factor(lab1, levels = 1:K1), median))
loc$iqr_x_km <- as.numeric(tapply(X[, 1], factor(lab1, levels = 1:K1), IQR))
loc$iqr_y_km <- as.numeric(tapply(X[, 2], factor(lab1, levels = 1:K1), IQR))
loc$main_station <- as.character(tapply(stead$covariates$station, factor(lab1, levels = 1:K1),
                                        function(s) names(which.max(table(s)))))
loc$share_main_station <- as.numeric(tapply(stead$covariates$station, factor(lab1, levels = 1:K1),
                                            function(s) max(table(s)) / length(s)))
write.csv(loc, "results/app_stead_type_locations.csv", row.names = FALSE)
print(loc, digits = 3)

# distance of the centre of the most concentrated type from Mount St. Helens
k_star <- which.max(types$max_Z)
cx <- median(X[lab1 == k_star, 1])
cy <- median(X[lab1 == k_star, 2])
d_msh <- sqrt((cx - volcanoes$x[1])^2 + (cy - volcanoes$y[1])^2)
cat("most concentrated type:", k_star, " centre", round(cx, 1), round(cy, 1),
    "km; distance from Mount St. Helens", round(d_msh, 1), "km\n")
# time span of each type (origin times, if available)
if (!is.null(m$source_origin_time)) {
  yrs <- as.numeric(substr(as.character(m$source_origin_time), 1, 4))
  print(tapply(yrs, factor(lab1, levels = 1:K1), range))
}

## ---- 2. map: 2 x 2 panels (types 1..K1, then the detail of the most concentrated type)
draw_panel <- function(k, xlim, ylim, letter, main, cex_pt) {
  plot(NA, xlim = xlim, ylim = ylim, asp = 1, xlab = "x (km)", ylab = "y (km)",
       bty = "n")
  panel_title(letter, main)
  points(X[, 1], X[, 2], pch = 16, cex = 0.25, col = bg_col)
  sel <- lab1 == k
  points(X[sel, 1], X[sel, 2], pch = type_pch[k], cex = cex_pt,
         col = transparent(type_cols[k], 0.75))
  points(stations_xy, pch = 25, bg = "white", cex = 0.7)
  points(volcanoes$x, volcanoes$y, pch = 24, bg = "red", cex = 1.0)
  box(lwd = lwd_box)
}
full_x <- range(X[, 1]); full_y <- range(X[, 2])
n_panels <- K1 + 1
open_pdf(file.path("..", "output", "figures", "fig_app_map.pdf"),
         height = ifelse(n_panels <= 4, 7.4, 10.5))
set_panels(c(ceiling(n_panels / 2), 2), mar = c(3, 3.3, 1.5, 0.4), oma = c(1.8, 0, 0, 0))
for (k in 1:K1) {
  draw_panel(k, full_x, full_y, letters[k],
             paste0("Type ", k, " (", sum(lab1 == k), " events)"), 0.4)
}
half <- 15
draw_panel(k_star, cx + c(-half, half), cy + c(-half, half), letters[K1 + 1],
           paste0("Type ", k_star, ", detail"), 0.7)
text(volcanoes$x, volcanoes$y, volcanoes$name, pos = 2, cex = 0.85, font = 2, xpd = NA)
legend_bottom(legend = c(paste("type", 1:K1), "all events", "station", "volcano"),
              pch = c(type_pch[1:K1], 16, 25, 24),
              col = c(type_cols[1:K1], bg_col, "black", "black"),
              pt.bg = c(rep(NA, K1 + 1), "white", "red"),
              pt.cex = c(rep(1, K1), 1, 0.8, 1), ncol = K1 + 3)
dev.off()

## ---- 3. spectra and example waveforms --------------------------------------------------
open_pdf(file.path("..", "output", "figures", "fig_app_curves.pdf"), height = 2.9)
set_panels(c(1, 2), mar = c(2.9, 3.3, 1.5, 0.5), oma = c(1.6, 0, 0, 0))
spec_max <- max(sapply(1:K1, function(k) max(colMeans(Z[lab1 == k, , drop = FALSE]))))
plot(NA, xlim = c(0.5, max(freq)), ylim = c(0, 1.05 * spec_max), log = "x", bty = "n",
     xlab = "frequency (Hz)", ylab = "normalised modulus")
panel_title("a", "Mean spectrum by type")
for (k in 1:K1) {
  lines(freq, colMeans(Z[lab1 == k, , drop = FALSE]), col = type_cols[k], lwd = lwd_main)
}
box(lwd = lwd_box)
tt <- (1:ncol(Y)) / stead$fs - stead$pre_p_sec
plot(NA, xlim = range(tt), ylim = c(-K1 + 0.5, 0.6), yaxt = "n", bty = "n",
     xlab = "time after the P arrival (s)", ylab = "")
panel_title("b", "Representative waveforms")
for (k in 1:K1) {
  idx <- which(lab1 == k)
  med <- apply(Z[idx, , drop = FALSE], 2, median)
  d <- rowSums((Z[idx, , drop = FALSE] - matrix(med, nrow = length(idx), ncol = ncol(Z), byrow = TRUE))^2)
  y <- Y[idx[which.min(d)], ]
  lines(tt, 0.45 * y / max(abs(y)) - (k - 1), col = type_cols[k], lwd = 0.7)
}
axis(2, at = -(0:(K1 - 1)), labels = paste("type", 1:K1), tick = FALSE)
box(lwd = lwd_box)
legend_bottom(legend = paste("type", 1:K1), col = type_cols[1:K1], lwd = lwd_main,
              ncol = K1)
dev.off()
cat("written fig_app_map.pdf and fig_app_curves.pdf\n")
