################################################################################
#  app_04_competitors_stead.R
#
#  Competing tests on the STEAD application (Section 7 of the paper), applied
#  to the same events, the same primary clustering and the same strata as
#  FCST (app_02_fcst_stead.R must be run first).
#
#    on the estimated types (primary clustering):
#      segregation test (Diggle et al. 2005), multinomial scan (Jung et al. 2010)
#    on the curves (Fourier-modulus representation, no clustering):
#      functional scan, mark-weighted K-function, distance covariance
#
#  Each test is run with unrestricted permutations and with permutations
#  within the distance x magnitude strata. For the multinomial scan the most
#  likely window is also saved (centre and number of events), to be compared
#  with the types found concentrated by FCST.
#
#  Output: results/app_stead_competitors.csv
#          ../output/tables/tab_app_competitors.tex
################################################################################

source("fcst_methods.R")
source("fcst_competitors.R")
set.seed(20261120)

B        <- 999        # the scan statistics are slow with n ~ 3500
max_frac <- 0.25       # largest scan window: a quarter of the events

stead <- readRDS("results/app_stead_data.rds")
lab   <- readRDS("results/app_stead_labels.rds")
X <- stead$X
Y <- stead$Y
labels <- lab$labels[, lab$primary]
Z <- represent_curves(Y, method = "spectrum")
radii <- lab$radii
n <- nrow(X)
cat("events:", n, " primary clustering:", lab$primary, "\n")

O <- neighbour_order(X, max_frac)            # shared by the two scans
pairs <- build_pair_list(X, radii, kernel = "ball")

run_all <- function(strata, label) {
  cat("\n== competing tests:", label, "==\n")
  out <- list()
  tm <- system.time(r <- fcst_test(X, labels, radii, B = B, strata = strata,
                                   progress = TRUE, progress_label = "1/6 FCST"))
  out$fcst <- c(r$p_global, tm["elapsed"])
  tm <- system.time(r <- segregation_test(X, labels, B = B, strata = strata,
                                          progress = TRUE, progress_label = "2/6 segregation"))
  out$segregation <- c(r$p_value, tm["elapsed"])
  tm <- system.time(r_scan <- multinomial_scan(X, labels, B = B, O = O, strata = strata,
                                               progress = TRUE, progress_label = "3/6 multinomial scan"))
  out$mscan <- c(r_scan$p_value, tm["elapsed"])
  tm <- system.time(r <- functional_scan(X, Z, B = B, O = O, strata = strata,
                                         progress = TRUE, progress_label = "4/6 functional scan"))
  out$fscan <- c(r$p_value, tm["elapsed"])
  tm <- system.time(r <- mark_weighted_K(X, Z, radii, B = B, pairs = pairs, strata = strata,
                                         progress = TRUE, progress_label = "5/6 mark-weighted K"))
  out$markK <- c(r$p_value, tm["elapsed"])
  tm <- system.time(r <- dcov_test(X, Z, B = B, strata = strata,
                                   progress = TRUE, progress_label = "6/6 distance covariance"))
  out$dcov <- c(r$p_value, tm["elapsed"])
  res <- data.frame(analysis = label, test = names(out),
                    p_value = sapply(out, function(v) v[1]),
                    seconds = sapply(out, function(v) v[2]))
  # most likely window of the multinomial scan
  members <- O[r_scan$window_centre, 1:r_scan$window_size]
  res$scan_window_n <- r_scan$window_size
  res$scan_window_x <- mean(X[members, 1])
  res$scan_window_y <- mean(X[members, 2])
  res$scan_window_type_share <- paste(round(table(factor(labels[members],
                                       levels = sort(unique(labels)))) /
                                       length(members), 2), collapse = "/")
  return(res)
}

res <- rbind(run_all(NULL, "unrestricted"),
             run_all(lab$strata_D, "distance x magnitude"))
write.csv(res, "results/app_stead_competitors.csv", row.names = FALSE)
print(res, digits = 3)

## ---- table for the paper -----------------------------------------------------------
fmt_p <- function(p) ifelse(p < 0.001, "$<0.001$", formatC(p, format = "f", digits = 3))
test_label <- c(fcst = "FCST", segregation = "Segregation", mscan = "Multinomial scan",
                fscan = "Functional scan", markK = "Mark-weighted $K$",
                dcov = "Distance covariance")
input_label <- c(fcst = "types", segregation = "types", mscan = "types",
                 fscan = "curves", markK = "curves", dcov = "curves")
tests <- unique(res$test)
lines <- c("\\begin{tabular}{@{}llrrr@{}}", "\\toprule",
           "Test & Input & Unrestricted & Distance $\\times$ magnitude & Time (s) \\\\",
           "\\midrule")
for (t in tests) {
  pu <- res$p_value[res$test == t & res$analysis == "unrestricted"]
  pd <- res$p_value[res$test == t & res$analysis == "distance x magnitude"]
  sec <- res$seconds[res$test == t & res$analysis == "unrestricted"]
  lines <- c(lines, paste0(test_label[[t]], " & ", input_label[[t]], " & ", fmt_p(pu),
                           " & ", fmt_p(pd), " & ", formatC(sec, format = "f", digits = 0),
                           " \\\\"))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path("..", "paper", "tables", "tab_app_competitors.tex"))
cat("written ../output/tables/tab_app_competitors.tex\n")
