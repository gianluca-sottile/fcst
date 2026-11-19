################################################################################
#  app_02_fcst_stead.R
#
#  Application (Section 7 of the paper): FCST on the STEAD selection prepared
#  by app_01_prepare_stead.R (results/app_stead_data.rds).
#
#  Analysis plan (fixed before looking at the results)
#    Representation : modulus of the Fourier transform of the vertical
#                     component, rescaled to unit norm (invariant to
#                     amplitude, hence to magnitude and geometrical spreading,
#                     and to the phase).
#    Clusterings    : k-means and Gaussian mixture on the FPCA scores
#                     (90% of variance), K = 2, 3, 4  ->  A = 6 clusterings,
#                     calibrated jointly (Section 3.4 of the paper).
#                     K = 2, ..., 5 (8 clusterings calibrated jointly).
#                     Primary clustering for maps and description: Gaussian
#                     mixture with K = 3, at the elbow of the ICL and BIC
#                     curves computed from the curves only (Assumption 1
#                     still holds); the criteria are printed for the record.
#    Radii          : 2, 5, 10, 20 km.
#    Tests          : (U) unrestricted;
#                     (D) strata = tertiles of source-receiver distance x
#                         tertiles of magnitude (9 strata);
#                     (S) strata = recording station (stations with fewer
#                         than 'min_station' events are pooled), only when
#                         trace_mode = "nearest".
#    B              : 1999 permutations (9999 for the final run).
#
#  Outputs
#    results/app_stead_tests.csv        global and type-specific results
#    results/app_stead_types.csv        description of the types
#    ../output/tables/tab_app_tests.tex  table for the paper
#    ../output/tables/tab_app_types.tex  table for the paper
#    ../output/figures/fig_app_map.pdf   map of the events by type
#    ../output/figures/fig_app_curves.pdf mean spectra and example waveforms
################################################################################

source("fcst_methods.R")
set.seed(20261119)

dir_tab <- file.path("..", "output", "tables")
dir_fig <- file.path("..", "output", "figures")

## ---- settings ----------------------------------------------------------------------
B           <- 1999                          # 9999 for the final run if time allows
radii       <- c(2, 5, 10, 20)               # km
K_values    <- c(2, 3, 4, 5)                  # same range for the joint family and the choice of K
K_primary   <- 3                             # see the ICL/BIC/silhouette summary printed below
algorithms  <- c("fpca_kmeans", "fpca_mclust")
n_tertiles  <- 3
min_station <- 10
alpha       <- 0.05

stead <- readRDS("results/app_stead_data.rds")
X  <- stead$X
Y  <- stead$Y
cv <- stead$covariates
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
for (v in c("distance_km", "magnitude", "depth_km")) cv[[v]] <- num(cv[[v]])
if (anyNA(cv$distance_km) || anyNA(cv$magnitude)) {
  stop("missing distance or magnitude: needed for the strata")
}
n  <- nrow(X)
cat("events:", n, "\n")

## ---- 1. clusterings ----------------------------------------------------------------
Z <- represent_curves(Y, method = "spectrum")
typing_names <- character(0)
labels <- NULL
for (alg in algorithms) {
  for (K in K_values) {
    cat("clustering:", alg, "K =", K, "\n")
    lab <- type_curves(Y, K, representation = "spectrum", algorithm = alg)
    labels <- cbind(labels, lab)
    typing_names <- c(typing_names, paste0(alg, "_K", K))
  }
}
colnames(labels) <- typing_names
# primary clustering: criteria computed from the FPCA scores of the curves
# only (no locations), printed for the record
S_scores <- fpca_scores(Z, 0.90)
bic_fit <- mclust::Mclust(S_scores, G = K_values, verbose = FALSE)
icl_fit <- mclust::mclustICL(S_scores, G = K_values, verbose = FALSE)
crit <- data.frame(K = K_values,
                   BIC = apply(bic_fit$BIC, 1, max, na.rm = TRUE),
                   ICL = apply(icl_fit, 1, max, na.rm = TRUE),
                   silhouette_kmeans = NA)
d_scores <- dist(S_scores)
for (j in seq_along(K_values)) {
  km <- kmeans(S_scores, K_values[j], nstart = 25)
  crit$silhouette_kmeans[j] <- mean(cluster::silhouette(km$cluster, d_scores)[, 3])
}
crit$ICL_gain <- c(NA, diff(crit$ICL))
print(crit, digits = 6)
write.csv(crit, "results/app_stead_choice_K.csv", row.names = FALSE)
primary <- list(algorithm = "fpca_mclust", K = K_primary)
cat("primary clustering: Gaussian mixture with K =", primary$K, "\n")
primary_name <- paste0(primary$algorithm, "_K", primary$K)
a_primary <- which(typing_names == primary_name)
cat("cluster sizes, primary clustering:", table(labels[, a_primary]), "\n")

## ---- 2. strata ---------------------------------------------------------------------
strata_D <- make_strata(cbind(cv$distance_km, cv$magnitude), n_bins = n_tertiles)
strata_S <- NULL
if (stead$trace_mode == "nearest") {
  st <- as.character(cv$station)
  small <- names(which(table(st) < min_station))
  st[st %in% small] <- "pooled"
  strata_S <- as.integer(factor(st))
}

## ---- 3. tests ----------------------------------------------------------------------
run_analysis <- function(strata, label) {
  cat("\n== test:", label, "(", ncol(labels), "clusterings,", B, "permutations ) ==\n")
  test <- fcst_test(X, labels, radii, B = B, strata = strata,
                    progress = TRUE, progress_label = label)
  tab <- test$table
  tab$clustering <- typing_names[tab$typing]
  tab$analysis <- label
  glob <- data.frame(analysis = label, clustering = typing_names,
                     p_clustering = test$p_typing, p_joint = test$p_global)
  return(list(types = tab, global = glob))
}
res <- list(U = run_analysis(NULL, "unrestricted"),
            D = run_analysis(strata_D, "distance x magnitude"))
if (!is.null(strata_S)) res$S <- run_analysis(strata_S, "station")

types_all  <- do.call(rbind, lapply(res, function(r) r$types))
global_all <- do.call(rbind, lapply(res, function(r) r$global))
write.csv(types_all,  "results/app_stead_tests_types.csv",  row.names = FALSE)
write.csv(global_all, "results/app_stead_tests_global.csv", row.names = FALSE)
saveRDS(list(labels = labels, typing_names = typing_names, primary = primary_name,
             K_primary = primary$K, criteria = crit,
             strata_D = strata_D, strata_S = strata_S, radii = radii),
        "results/app_stead_labels.rds")   # used by app_04_competitors_stead.R
print(global_all, digits = 3)

## ---- 4. scale at which each type is most concentrated (primary clustering) ----------
pairs <- build_pair_list(X, radii, kernel = "ball")
K1 <- primary$K
lab1 <- labels[, a_primary]
mom <- all_null_moments(lab1, pairs, K1)
Zobs <- standardise_pair(pair_statistic(lab1, pairs, K1), mom, K1)
best_radius <- radii[apply(Zobs, 1, which.max)]

## ---- 5. description of the types (primary clustering) -------------------------------
freq <- (1:ncol(Z)) * stead$fs / ncol(Y)          # Hz
describe <- data.frame(type = 1:K1, n = as.integer(table(factor(lab1, levels = 1:K1))))
describe$peak_freq_hz <- NA
describe$median_distance_km <- NA
describe$median_magnitude <- NA
describe$median_depth_km <- NA
for (k in 1:K1) {
  sel <- lab1 == k
  mean_spec <- colMeans(Z[sel, , drop = FALSE])
  describe$peak_freq_hz[k] <- freq[which.max(mean_spec)]
  describe$median_distance_km[k] <- median(cv$distance_km[sel])
  describe$median_magnitude[k] <- median(cv$magnitude[sel])
  describe$median_depth_km[k] <- median(cv$depth_km[sel], na.rm = TRUE)
}
describe$best_radius_km <- best_radius
describe$max_Z <- apply(Zobs, 1, max)
for (an in names(res)) {
  tt <- res[[an]]$types
  tt <- tt[tt$clustering == primary_name, ]
  describe[[paste0("p_holm_", an)]] <- tt$p_holm[match(1:K1, tt$type)]
}
write.csv(describe, "results/app_stead_types.csv", row.names = FALSE)
print(describe, digits = 3)

## ---- 6. tables for the paper ---------------------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))
fmt_p <- function(p) ifelse(p < 0.001, "$<0.001$", fmt(p))
alg_label <- c(fpca_kmeans = "$k$-means", fpca_mclust = "Gaussian mixture")

lines <- c(paste0("\\begin{tabular}{@{}ll", strrep("r", length(res)), "@{}}"),
           "\\toprule",
           paste0("Clustering & $K$ & ",
                  paste(c(U = "Unrestricted", D = "Distance $\\times$ magnitude",
                          S = "Station")[names(res)], collapse = " & "), " \\\\"),
           "\\midrule")
for (a in seq_along(typing_names)) {
  alg <- sub("_K[0-9]+$", "", typing_names[a])
  K <- sub("^.*_K", "", typing_names[a])
  vals <- sapply(res, function(r) r$global$p_clustering[a])
  lines <- c(lines, paste0(alg_label[[alg]], " & ", K, " & ",
                           paste(fmt_p(vals), collapse = " & "), " \\\\"))
}
joint <- sapply(res, function(r) r$global$p_joint[1])
lines <- c(lines, "\\midrule",
           paste0("\\multicolumn{2}{@{}l}{Joint test} & ",
                  paste(fmt_p(joint), collapse = " & "), " \\\\"),
           "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(dir_tab, "tab_app_tests.tex"))

lines <- c(paste0("\\begin{tabular}{@{}rrrrrrr", strrep("r", length(res)), "@{}}"),
           "\\toprule",
           paste0("Type & $n_k$ & Peak (Hz) & Distance (km) & Magnitude & Depth (km) & ",
                  "Radius (km) & ",
                  paste(paste0("$p$ (", names(res), ")"), collapse = " & "), " \\\\"),
           "\\midrule")
for (k in 1:K1) {
  pv <- sapply(names(res), function(an) describe[[paste0("p_holm_", an)]][k])
  lines <- c(lines, paste0(k, " & ", describe$n[k], " & ",
                           fmt(describe$peak_freq_hz[k], 1), " & ",
                           fmt(describe$median_distance_km[k], 1), " & ",
                           fmt(describe$median_magnitude[k], 2), " & ",
                           fmt(describe$median_depth_km[k], 1), " & ",
                           describe$best_radius_km[k], " & ",
                           paste(fmt_p(pv), collapse = " & "), " \\\\"))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(dir_tab, "tab_app_types.tex"))

## ---- 7. figures --------------------------------------------------------------------
type_cols <- c("#D55E00", "#0072B2", "#7F7F7F", "#009E73", "#CC79A7")
type_pch  <- c(16, 17, 15, 18, 3)

# volcanoes in the window (for orientation only)
volcanoes <- data.frame(
  name = c("Mount St. Helens", "Mount Hood", "Mount Adams", "Mount Rainier"),
  lat  = c(46.1912, 45.3736, 46.2024, 46.8528),
  lon  = c(-122.1944, -121.6960, -121.4909, -121.7604))
vx <- (volcanoes$lon - stead$lon0) * stead$km_per_deg["lon"]
vy <- (volcanoes$lat - stead$lat0) * stead$km_per_deg["lat"]

pdf(file.path(dir_fig, "fig_app_map.pdf"), width = 10, height = 5.2)
par(mfrow = c(1, K1), mar = c(3.5, 3.5, 2.5, 0.5), mgp = c(2.2, 0.7, 0))
for (k in 1:K1) {
  plot(X[, 1], X[, 2], pch = 16, cex = 0.35, col = "grey80", asp = 1,
       xlab = "x (km)", ylab = "y (km)",
       main = paste0("Type ", k, " (n = ", sum(lab1 == k), ")"))
  sel <- lab1 == k
  points(X[sel, 1], X[sel, 2], pch = type_pch[k], col = type_cols[k], cex = 0.6)
  points(vx, vy, pch = 24, bg = "white", cex = 1.3)
  points(unique(stead$station_xy), pch = 2, cex = 0.6, col = "black")
}
dev.off()

pdf(file.path(dir_fig, "fig_app_curves.pdf"), width = 10, height = 4.2)
par(mfrow = c(1, 2), mar = c(3.8, 3.8, 2.5, 0.8), mgp = c(2.3, 0.7, 0))
plot(NA, xlim = range(freq), ylim = c(0, max(Z)), log = "x",
     xlab = "frequency (Hz)", ylab = "normalised Fourier modulus",
     main = "Mean spectrum by type")
for (k in 1:K1) {
  lines(freq, colMeans(Z[lab1 == k, , drop = FALSE]), col = type_cols[k], lwd = 2.5)
}
legend("topright", bty = "n", legend = paste("type", 1:K1),
       col = type_cols[1:K1], lwd = 2.5)
tt <- (1:ncol(Y)) / stead$fs - stead$pre_p_sec
plot(NA, xlim = range(tt), ylim = c(-K1 + 0.5, 0.6), yaxt = "n",
     xlab = "time after P arrival (s)", ylab = "",
     main = "Example waveforms (normalised)")
for (k in 1:K1) {
  idx <- which(lab1 == k)
  # the event closest to the type-wise median spectrum
  d <- rowSums((Z[idx, , drop = FALSE] -
                matrix(apply(Z[idx, , drop = FALSE], 2, median),
                       nrow = length(idx), ncol = ncol(Z), byrow = TRUE))^2)
  y <- Y[idx[which.min(d)], ]
  lines(tt, 0.45 * y / max(abs(y)) - (k - 1), col = type_cols[k])
}
axis(2, at = -(0:(K1 - 1)), labels = paste("type", 1:K1), las = 1, tick = FALSE)
dev.off()
cat("tables and figures written\n")
