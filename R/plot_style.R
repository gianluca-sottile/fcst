################################################################################
#  plot_style.R
#
#  Common graphical settings for all the figures of the paper
#  (make_tables_figures.R, make_design_figures.R, app_03_data_map.R,
#  app_05_figures_stead.R).
#
#  Principle: every figure is drawn at its FINAL printed size (the text width
#  of the paper is about 6.5 inches), with a 9-point font, so that what is
#  drawn in R is what is seen in the PDF of the paper, without reductions.
################################################################################

## ---- sizes ---------------------------------------------------------------------------
fig_width   <- 6.5      # inches, full text width
font_size   <- 9        # points

## ---- colours (Okabe-Ito palette: distinguishable also by colour-blind readers) -------
col_types <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00")  # shape types
pch_types <- c(16, 17, 15, 18, 8)
col_grey_bg   <- "grey82"   # background points (all events)
col_reference <- "grey45"   # nominal level, contours, reference lines

# accuracy of the clustering (power figure): dark to light, all clearly visible
col_accuracy <- c("#1B1B1B", "#0072B2", "#009E73", "#CC79A7")
pch_accuracy <- c(16, 17, 15, 18)

# tests (comparison with competitors): solid = tests on types,
# dashed = tests on curves; FCST black and thicker
tests_all  <- c("fcst", "centroid", "segregation", "mscan", "fscan", "markK", "dcov")
tests_name <- c(fcst = "FCST", centroid = "centroid", segregation = "segregation",
                mscan = "multinomial scan", fscan = "functional scan",
                markK = "mark-weighted K", dcov = "distance covariance")
col_tests <- c(fcst = "#000000", centroid = "#999999", segregation = "#D55E00",
               mscan = "#E69F00", fscan = "#0072B2", markK = "#56B4E9",
               dcov = "#009E73")
lty_tests <- c(fcst = 1, centroid = 3, segregation = 1, mscan = 1, fscan = 2,
               markK = 2, dcov = 2)
pch_tests <- c(fcst = 16, centroid = 4, segregation = 17, mscan = 15, fscan = 1,
               markK = 2, dcov = 0)
lwd_tests <- c(fcst = 2.2, centroid = 1.4, segregation = 1.4, mscan = 1.4,
               fscan = 1.4, markK = 1.4, dcov = 1.4)

## ---- line widths ---------------------------------------------------------------------
lwd_main <- 1.6
lwd_thin <- 0.6
lwd_box  <- 0.8

## ---- helpers -------------------------------------------------------------------------
open_pdf <- function(file, height, width = fig_width) {
  pdf(file, width = width, height = height, pointsize = font_size,
      family = "Helvetica", useDingbats = FALSE)
}

set_panels <- function(mfrow, mar = c(3, 3.2, 1.6, 0.5), oma = c(0, 0, 0, 0)) {
  # NB: par(mfrow) shrinks the text when there are 3 or more rows or columns;
  # cex = 1 is restored afterwards so that every text is at font_size points.
  par(mfrow = mfrow)
  par(cex = 1, mar = mar, oma = oma, mgp = c(1.8, 0.45, 0), tcl = -0.25,
      cex.main = 1, font.main = 1, cex.lab = 1, cex.axis = 0.9, las = 1,
      lwd = lwd_box)
}

panel_title <- function(letter, text) {
  # "(a) text" on the left of the panel title line
  title(main = paste0("(", letter, ") ", text), adj = 0, line = 0.5)
}

legend_bottom <- function(..., ncol = 3) {
  # legend in the bottom outer margin, across the whole figure
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("bottom", bty = "n", xpd = NA, ncol = ncol, cex = 1,
         inset = c(0, 0.005), x.intersp = 0.6, seg.len = 1.8, ...)
}

transparent <- function(col, alpha) adjustcolor(col, alpha.f = alpha)
