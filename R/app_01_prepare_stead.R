################################################################################
#  app_01_prepare_stead.R
#
#  Step 1 of the STEAD application: combines the per-chunk extracts written by
#  app_00_extract_chunk.R (results/stead_extract_chunkN.rds), keeps ONE trace
#  per event and builds the data used by app_02_fcst_stead.R.
#  The original STEAD files are no longer needed at this stage.
#
#  Selection (applied in app_00): local earthquakes, 45-47 N, 123-121 W,
#  receiver_type "EH", source_magnitude_type "MD"; vertical component, window
#  from 0.5 s before to 12 s after the P arrival, detrended, tapered and
#  band-pass filtered (1-20 Hz).
#
#  One trace per event ('trace_mode'):
#    "nearest" : the trace recorded at the station closest to the epicentre
#    "station" : only events recorded at one station (the one with most
#                events, or 'station_code'): same receiver for all curves
#
#  Output: results/app_stead_data.rds, results/app_stead_station_counts.csv
################################################################################

library(data.table)

## ---- settings ----------------------------------------------------------------------
# extracts to combine: the SeisBench extract if present (app_00_extract_seisbench.R),
# otherwise the per-chunk extracts (app_00_extract_chunk.R)
if (file.exists(file.path("results", "stead_extract_seisbench.rds"))) {
  chunks <- "seisbench"
} else {
  chunks <- c("chunk2", "chunk3", "chunk4", "chunk5", "chunk6")
}
trace_mode   <- "nearest"
station_code <- NULL
lat_range    <- c(45, 47)
lon_range    <- c(-123, -121)

## ---- 1. combine the extracts -----------------------------------------------------------
meta_list <- list()
Y_list <- list()
for (ch in chunks) {
  f <- file.path("results", paste0("stead_extract_", ch, ".rds"))
  if (!file.exists(f)) {
    warning("missing extract: ", f, " (run app_00_extract_chunk.R ", ch, ")")
    next
  }
  e <- readRDS(f)
  meta_list[[ch]] <- e$meta
  Y_list[[ch]] <- e$Y
  settings <- e[c("fs", "pre_p_sec", "post_p_sec", "band_hz", "filtered")]
}
meta <- rbindlist(meta_list, fill = TRUE)
Y_all <- do.call(rbind, Y_list)
meta[, row_id := .I]
cat("traces:", nrow(meta), " events:", uniqueN(meta$source_id),
    " stations:", uniqueN(meta$receiver_code), "\n")

station_counts <- meta[, .(n_events = uniqueN(source_id),
                           lat = receiver_latitude[1],
                           lon = receiver_longitude[1]),
                       by = .(network_code, receiver_code)][order(-n_events)]
fwrite(station_counts, "results/app_stead_station_counts.csv")
print(head(station_counts, 10))

## ---- 2. one trace per event ------------------------------------------------------------
if (trace_mode == "nearest") {
  setorder(meta, source_id, source_distance_km)
  one <- meta[, .SD[1], by = source_id]
} else if (trace_mode == "station") {
  if (is.null(station_code)) station_code <- station_counts$receiver_code[1]
  one <- meta[receiver_code == station_code]
  setorder(one, source_id, source_distance_km)
  one <- one[, .SD[1], by = source_id]
  cat("station:", station_code, "\n")
} else {
  stop("unknown trace_mode")
}
Y <- Y_all[one$row_id, , drop = FALSE]
cat("events kept:", nrow(one), "\n")

## ---- 3. coordinates (km) and covariates ------------------------------------------------
lat0 <- mean(lat_range)
lon0 <- mean(lon_range)
km_per_deg_lat <- 110.57
km_per_deg_lon <- 111.32 * cos(lat0 * pi / 180)
X <- cbind(x = (one$source_longitude - lon0) * km_per_deg_lon,
           y = (one$source_latitude  - lat0) * km_per_deg_lat)
station_xy <- cbind(x = (one$receiver_longitude - lon0) * km_per_deg_lon,
                    y = (one$receiver_latitude  - lat0) * km_per_deg_lat)
# some metadata columns are read as text (e.g. missing values written as
# "None"): convert explicitly; non-numeric entries become NA
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
covariates <- data.frame(distance_km = num(one$source_distance_km),
                         magnitude   = num(one$source_magnitude),
                         depth_km    = num(one$source_depth_km),
                         station     = one$receiver_code,
                         network     = one$network_code)

stead <- c(list(X = X, Y = Y, covariates = covariates, meta = one,
                station_xy = station_xy, trace_mode = trace_mode,
                lat0 = lat0, lon0 = lon0,
                km_per_deg = c(lon = km_per_deg_lon, lat = km_per_deg_lat),
                n_traces = nrow(meta), n_stations = uniqueN(meta$receiver_code)),
           settings)
saveRDS(stead, "results/app_stead_data.rds")
cat("missing values: distance", sum(is.na(covariates$distance_km)),
    " magnitude", sum(is.na(covariates$magnitude)),
    " depth", sum(is.na(covariates$depth_km)), "\n")
cat("saved results/app_stead_data.rds:", nrow(Y), "events x", ncol(Y), "samples\n")
