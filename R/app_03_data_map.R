################################################################################
#  app_03_data_map.R
#
#  Figure describing the STEAD selection (Section 7.1 of the paper):
#    (a) epicentres in the whole window, coloured by depth, with stations and
#        volcanoes (state borders if the package 'maps' is installed);
#    (b) zoom on Mount St. Helens, where most events occur;
#    (c) magnitude against source-receiver distance.
#  Output: ../output/figures/fig_app_data.pdf
#  Input : results/app_stead_data.rds (app_01_prepare_stead.R)
################################################################################

source("plot_style.R")
stead <- readRDS("results/app_stead_data.rds")
m  <- stead$meta
cv <- stead$covariates
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
for (v in c("distance_km", "magnitude", "depth_km")) cv[[v]] <- num(cv[[v]])
have_maps <- requireNamespace("maps", quietly = TRUE)

volcanoes <- data.frame(
  name = c("St. Helens", "Hood", "Adams", "Rainier"),
  lat  = c(46.1912, 45.3736, 46.2024, 46.8528),
  lon  = c(-122.1944, -121.6960, -121.4909, -121.7604))
stations <- unique(data.frame(lon = m$receiver_longitude, lat = m$receiver_latitude))

# colour by depth (km), 6 classes
depth_breaks <- c(-Inf, 2, 4, 6, 10, 20, Inf)
depth_labels <- c("< 2", "2-4", "4-6", "6-10", "10-20", "> 20")
depth_cols   <- c("#E8C200", "#7AD151", "#22A884", "#2A788E", "#414487", "#440154")
depth_class  <- cut(cv$depth_km, depth_breaks, labels = FALSE)
pt_col <- adjustcolor(depth_cols[depth_class], alpha.f = 0.7)
pt_col[is.na(depth_class)] <- "grey50"          # depth not available

asp_lonlat <- 1 / cos(46 * pi / 180)      # 1 degree of longitude is shorter

draw_map <- function(xlim, ylim, letter, main, cex_pt) {
  plot(NA, xlim = xlim, ylim = ylim, asp = asp_lonlat, xlab = "longitude",
       ylab = "latitude", bty = "n")
  panel_title(letter, main)
  if (have_maps) maps::map("state", add = TRUE, col = "grey55", lwd = 0.7)
  points(m$source_longitude, m$source_latitude, pch = 16, cex = cex_pt, col = pt_col)
  points(stations$lon, stations$lat, pch = 25, bg = "white", col = "black", cex = 0.8)
  points(volcanoes$lon, volcanoes$lat, pch = 24, bg = "red", col = "black", cex = 1.1)
  inside <- volcanoes$lon >= xlim[1] & volcanoes$lon <= xlim[2] &
            volcanoes$lat >= ylim[1] & volcanoes$lat <= ylim[2]
  text(volcanoes$lon[inside], volcanoes$lat[inside], volcanoes$name[inside],
       pos = 4, cex = 0.85, font = 2, xpd = NA)
  box(lwd = lwd_box)
}

# layout: the regional map on the left (two rows), the detail and the
# magnitude-distance plot on the right
open_pdf(file.path("..", "paper", "figures", "fig_app_data.pdf"), height = 5.6)
layout(matrix(c(1, 2, 1, 3), nrow = 2, byrow = TRUE), widths = c(1.15, 1))
par(cex = 1, mar = c(3, 3.3, 1.5, 0.5), oma = c(2.4, 0, 0, 0), mgp = c(1.9, 0.45, 0),
    tcl = -0.25, cex.main = 1, font.main = 1, cex.axis = 0.85, las = 1)
draw_map(c(-123, -121), c(45, 47), "a", "Epicentres and stations", 0.35)
rect(-122.35, 46.10, -122.05, 46.30, border = "black", lty = 2)
draw_map(c(-122.35, -122.05), c(46.10, 46.30), "b", "Mount St. Helens", 0.5)
plot(cv$distance_km, cv$magnitude, pch = 16, cex = 0.4, col = pt_col, log = "x",
     xlab = "source-receiver distance (km)", ylab = "duration magnitude", bty = "n")
panel_title("c", "Magnitude and distance")
box(lwd = lwd_box)
legend_bottom(legend = c(paste(depth_labels, "km"), "station", "volcano"),
              pch = c(rep(16, 6), 25, 24), col = c(depth_cols, "black", "black"),
              pt.bg = c(rep(NA, 6), "white", "red"), pt.cex = c(rep(1.2, 6), 0.9, 1.1),
              ncol = 4)
dev.off()
cat("written ../output/figures/fig_app_data.pdf\n")

# numbers quoted in the paper
cat("events:", nrow(m), "\n")
cat("within 15 km of Mount St. Helens:",
    sum(sqrt(((m$source_longitude - volcanoes$lon[1]) * cos(46 * pi / 180))^2 +
             (m$source_latitude - volcanoes$lat[1])^2) * 111.2 < 15), "\n")
print(summary(cv[, c("distance_km", "magnitude", "depth_km")]))
print(table(cv$station))
