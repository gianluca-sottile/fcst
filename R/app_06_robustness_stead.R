################################################################################
#  app_06_robustness_stead.R
#
#  Robustness of the STEAD results (Section 7.3 and Section S9 of the
#  Supplementary Material). Requires app_00 (extract), app_01 and app_02.
#
#  Does the concentration of types 2 and 3 of the primary clustering survive
#    (a) strata that also include TIME: station x year, and distance x magnitude
#        x period (tertiles of each)? Type 3 was recorded by one station in
#        2013-2015, so time (state of the instrument, swarm) is the natural
#        alternative explanation;
#    (b) a SINGLE-STATION analysis (B201, the station with most events): all
#        curves share site and instrument;
#    (c) a shorter window (0.5 s before to 6 s after P), a log-spectrum
#        representation, finer strata (quintiles of distance x tertiles of
#        magnitude) and other radii?
#  (d) Waveform similarity: maximum normalised cross-correlation between pairs
#      of type-3 events and between pairs of type-1 events recorded by the same
#      station (B201): repeating earthquakes have correlations close to 1.
#
#  For analyses that re-estimate the clustering, each new type is matched to
#  primary types 2 and 3 by the Jaccard index of their sets of events.
#
#  Output: results/app_stead_robustness.csv, results/app_stead_xcorr.csv,
#          ../output/tables/tab_app_robustness.tex
################################################################################

source("fcst_methods.R")
set.seed(20261120)

B     <- 1999
radii <- c(2, 5, 10, 20)
K     <- 3
n_tertiles <- 3

stead <- readRDS("results/app_stead_data.rds")
lab   <- readRDS("results/app_stead_labels.rds")
X  <- stead$X
Y  <- stead$Y
m  <- as.data.frame(stead$meta)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
cv <- stead$covariates
for (v in c("distance_km", "magnitude", "depth_km")) cv[[v]] <- num(cv[[v]])
year <- num(substr(as.character(m$source_origin_time), 1, 4))
year[is.na(year)] <- median(year, na.rm = TRUE)     # (none expected)
primary <- lab$labels[, lab$primary]
ids <- as.character(m$source_id)
cat("events:", nrow(X), " years:", min(year, na.rm = TRUE), "-", max(year, na.rm = TRUE), "\n")

## ---- helpers -------------------------------------------------------------------------
analyse <- function(X, labels, strata, radii, label) {
  # FCST on one clustering; returns one row per type with Z_max, radius, Holm p
  cat("\n==", label, "==\n")
  test <- fcst_test(X, labels, radii, B = B, strata = strata,
                    progress = TRUE, progress_label = label)
  Kl <- max(labels)
  pairs <- build_pair_list(X, radii, kernel = "ball")
  Z <- standardise_pair(pair_statistic(labels, pairs, Kl),
                        all_null_moments(labels, pairs, Kl), Kl)
  data.frame(analysis = label, p_global = test$p_global, type = 1:Kl,
             n_k = as.integer(table(factor(labels, levels = 1:Kl))),
             Z_max = apply(Z, 1, max, na.rm = TRUE),
             radius = radii[apply(Z, 1, which.max)],
             p_holm = test$table$p_holm)
}

jaccard_match <- function(new_labels, new_ids, ref_type) {
  # type of new_labels best matching primary type 'ref_type' (on shared events)
  ref <- ids[primary == ref_type]
  best <- c(type = NA, jaccard = 0)
  for (k in sort(unique(new_labels[new_labels > 0]))) {
    a <- new_ids[new_labels == k]
    jac <- length(intersect(a, ref)) / length(union(a, intersect(ref, new_ids)))
    if (jac > best["jaccard"]) best <- c(type = k, jaccard = jac)
  }
  return(best)
}

add_matching <- function(res, new_labels, new_ids) {
  res$matches_primary <- NA
  res$jaccard <- NA
  for (t in c(2, 3)) {
    mt <- jaccard_match(new_labels, new_ids, t)
    if (!is.na(mt["type"])) {
      res$matches_primary[res$type == mt["type"]] <- t
      res$jaccard[res$type == mt["type"]] <- mt["jaccard"]
    }
  }
  return(res)
}

cluster_mixture <- function(Zrep) {
  S <- fpca_scores(Zrep, 0.90)
  relabel_canonical(cluster_mclust(S, K), Zrep)
}

log_spectrum <- function(Y) {
  Zs <- represent_curves(Y, method = "spectrum")
  L <- log(Zs + 1e-3)
  L <- L - rowMeans(L)
  L / sqrt(rowSums(L^2))
}

results <- list()

## ---- (a) strata including time ---------------------------------------------------------
st_year <- as.integer(factor(paste(cv$station, year)))
results$station_year <- analyse(X, primary, st_year, radii, "station x year")
dmt <- make_strata(cbind(cv$distance_km, cv$magnitude, year), n_bins = n_tertiles)
results$dist_mag_time <- analyse(X, primary, dmt, radii, "distance x magnitude x period")
for (r in c("station_year", "dist_mag_time")) {
  results[[r]]$matches_primary <- ifelse(results[[r]]$type %in% c(2, 3), results[[r]]$type, NA)
  results[[r]]$jaccard <- ifelse(results[[r]]$type %in% c(2, 3), 1, NA)
}

## ---- (c) finer strata and other radii (primary clustering) ----------------------------
fine <- make_strata(cbind(cv$distance_km), n_bins = 5) * 10 +
        make_strata(cbind(cv$magnitude), n_bins = 3)
results$fine_strata <- analyse(X, primary, fine, radii, "finer strata")
results$radii_small <- analyse(X, primary, NULL, c(1, 2, 5, 10), "radii 1-10 km")
results$radii_large <- analyse(X, primary, NULL, c(5, 10, 20, 40), "radii 5-40 km")
for (r in c("fine_strata", "radii_small", "radii_large")) {
  results[[r]]$matches_primary <- ifelse(results[[r]]$type %in% c(2, 3), results[[r]]$type, NA)
  results[[r]]$jaccard <- ifelse(results[[r]]$type %in% c(2, 3), 1, NA)
}

## ---- (c) re-estimated clusterings: shorter window, log spectrum -----------------------
n_short <- round((stead$pre_p_sec + 6) * stead$fs)
lab_short <- cluster_mixture(represent_curves(Y[, 1:n_short], method = "spectrum"))
results$short_window <- add_matching(
  analyse(X, lab_short, NULL, radii, "window -0.5 to 6 s"), lab_short, ids)
lab_log <- cluster_mixture(log_spectrum(Y))
results$log_spectrum <- add_matching(
  analyse(X, lab_log, NULL, radii, "log spectrum"), lab_log, ids)

## ---- (b) single station B201 ----------------------------------------------------------
ex <- readRDS("results/stead_extract_seisbench.rds")
sel <- which(as.character(ex$meta$receiver_code) == "B201")
sel <- sel[!duplicated(as.character(ex$meta$source_id[sel]))]
mb <- as.data.frame(ex$meta)[sel, ]
Yb <- ex$Y[sel, , drop = FALSE]
Xb <- cbind((num(mb$source_longitude) - stead$lon0) * stead$km_per_deg["lon"],
            (num(mb$source_latitude)  - stead$lat0) * stead$km_per_deg["lat"])
lab_b <- cluster_mixture(represent_curves(Yb, method = "spectrum"))
ids_b <- as.character(mb$source_id)
cat("single station B201:", nrow(Xb), "events\n")
results$station_B201 <- add_matching(
  analyse(Xb, lab_b, NULL, radii, "single station B201"), lab_b, ids_b)
strata_b <- make_strata(cbind(num(mb$source_distance_km), num(mb$source_magnitude)), n_bins = 3)
results$station_B201_strata <- add_matching(
  analyse(Xb, lab_b, strata_b, radii, "B201, distance x magnitude"), lab_b, ids_b)

all_res <- do.call(rbind, results)
write.csv(all_res, "results/app_stead_robustness.csv", row.names = FALSE)
print(all_res, digits = 3)

## ---- (d) waveform cross-correlation -------------------------------------------------
max_xcorr <- function(a, b, max_lag = 100) {
  max(abs(ccf(a, b, lag.max = max_lag, plot = FALSE)$acf))
}
pair_sample <- function(idx, n_pairs = 1500) {
  if (length(idx) < 2) return(matrix(integer(0), ncol = 2))
  p <- t(replicate(n_pairs, sample(idx, 2)))
  return(p)
}
at_b201 <- cv$station == "B201"
xc <- list()
for (k in 1:K) {
  idx <- which(primary == k & at_b201)
  pr <- pair_sample(idx)
  if (nrow(pr) > 0) {
    v <- apply(pr, 1, function(q) max_xcorr(Y[q[1], ], Y[q[2], ]))
    xc[[k]] <- data.frame(type = k, n_events_B201 = length(idx), median_xcorr = median(v),
                          share_above_0.7 = mean(v > 0.7), share_above_0.9 = mean(v > 0.9))
  }
}
xc <- do.call(rbind, xc)
write.csv(xc, "results/app_stead_xcorr.csv", row.names = FALSE)
print(xc, digits = 3)

## ---- table for the Supplementary Material ----------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))
fmt_p <- function(p) ifelse(is.na(p), "--", ifelse(p < 0.001, "$<0.001$", fmt(p, 4)))
nice <- c(station_year = "Strata: station $\\times$ year",
          dist_mag_time = "Strata: distance $\\times$ magnitude $\\times$ period",
          fine_strata = "Strata: distance (5) $\\times$ magnitude (3)",
          radii_small = "Radii 1--10 km", radii_large = "Radii 5--40 km",
          short_window = "Window $-0.5$ to $6$ s", log_spectrum = "Log spectrum",
          station_B201 = "Station B201",
          station_B201_strata = "Station B201, distance $\\times$ magnitude")
cells_for <- function(res, t) {
  r <- res[!is.na(res$matches_primary) & res$matches_primary == t, ]
  if (nrow(r) == 0) return("-- & -- & -- & --")
  jac <- if (r$jaccard[1] == 1) "--" else fmt(r$jaccard[1], 2)
  paste(fmt_p(r$p_holm[1]), fmt(r$Z_max[1], 1), r$radius[1], jac, sep = " & ")
}
lines <- c("\\begin{tabular}{@{}lrrrrrrrrr@{}}", "\\toprule",
           " & & \\multicolumn{4}{c}{Type 2 (low frequency)} & \\multicolumn{4}{c}{Type 3 (compact group)} \\\\",
           "\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}",
           "Analysis & $n$ & $p$ & $Z_{\\max}$ & $r$ & Jaccard & $p$ & $Z_{\\max}$ & $r$ & Jaccard \\\\",
           "\\midrule")
for (r in names(results)) {
  res <- results[[r]]
  lines <- c(lines, paste0(nice[[r]], " & ", sum(res$n_k), " & ", cells_for(res, 2), " & ",
                           cells_for(res, 3), " \\\\"))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path("..", "output", "tables", "tab_app_robustness.tex"))
cat("written tab_app_robustness.tex\n")

lines <- c("\\begin{tabular}{@{}lrrrr@{}}", "\\toprule",
           "Type & Events at B201 & Median & $>0.7$ & $>0.9$ \\\\", "\\midrule")
for (i in seq_len(nrow(xc))) {
  lines <- c(lines, paste0(xc$type[i], " & ", xc$n_events_B201[i], " & ",
                           fmt(xc$median_xcorr[i], 2), " & ", fmt(xc$share_above_0.7[i], 2),
                           " & ", fmt(xc$share_above_0.9[i], 2), " \\\\"))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path("..", "output", "tables", "tab_app_xcorr.tex"))
cat("written ../output/tables/tab_app_xcorr.tex\n")
