################################################################################
#  app_00_extract_seisbench.R
#
#  Step 0 of the STEAD application when the data were downloaded through
#  SeisBench (python: seisbench.data.STEAD()) instead of the original chunks.
#  It replaces app_00_extract_chunk.R and writes the same kind of output, so
#  that app_01 and app_02 are unchanged.
#
#  SeisBench layout (~/.seisbench/datasets/stead/)
#    metadata.csv    one row per trace; column names differ from the original
#                    STEAD names (they are renamed below)
#    waveforms.hdf5  traces stored in "buckets": trace_name "bucket17$0,:3,:6000"
#                    means row 0 (counting from 0) of dataset /data/bucket17,
#                    3 components x 6000 samples; the order of the components
#                    is stored in /data_format/component_order
#
#  Output: results/stead_extract_seisbench.rds
#  Requirements: data.table, rhdf5, signal (optional)
################################################################################

library(data.table)
library(rhdf5)
h5closeAll()   # close handles left open by a previous interrupted run
have_signal <- requireNamespace("signal", quietly = TRUE)

## ---- settings ----------------------------------------------------------------------
sb_dir     <- "~/.seisbench/datasets/stead"
lat_range  <- c(45, 47)
lon_range  <- c(-123, -121)
fs         <- 100      # Hz (checked against the file below)
pre_p_sec  <- 0.5
post_p_sec <- 12
band_hz    <- c(1, 20)
taper_frac <- 0.05

dir.create("results", showWarnings = FALSE)
meta_file <- path.expand(file.path(sb_dir, "metadata.csv"))
h5_file   <- path.expand(file.path(sb_dir, "waveforms.hdf5"))

## ---- 1. metadata: read, rename to the original STEAD names, select ------------------
meta <- fread(meta_file)
setnames(meta,
         old = c("station_network_code", "station_code", "trace_channel",
                 "station_latitude_deg", "station_longitude_deg",
                 "trace_p_arrival_sample", "trace_s_arrival_sample",
                 "source_latitude_deg", "source_longitude_deg"),
         new = c("network_code", "receiver_code", "receiver_type",
                 "receiver_latitude", "receiver_longitude",
                 "p_arrival_sample", "s_arrival_sample",
                 "source_latitude", "source_longitude"))
# Step-by-step selection, printing how many traces survive each condition
# (comparisons on text columns ignore case and surrounding spaces).
txt <- function(x) tolower(trimws(as.character(x)))
num <- function(x) suppressWarnings(as.numeric(x))
cat("all traces:", nrow(meta), "\n")
cat("values of trace_category:", paste(head(unique(meta$trace_category), 10), collapse = " | "), "\n")
sel <- meta[txt(trace_category) == "earthquake_local"]
cat("  earthquake_local:", nrow(sel), "\n")
sel <- sel[num(source_latitude) >= lat_range[1] & num(source_latitude) <= lat_range[2] &
           num(source_longitude) >= lon_range[1] & num(source_longitude) <= lon_range[2]]
cat("  in the window:", nrow(sel), "\n")
cat("  channels in the window:", paste(names(table(sel$receiver_type)), table(sel$receiver_type), collapse = ", "), "\n")
cat("  magnitude types in the window:", paste(names(table(sel$source_magnitude_type)), table(sel$source_magnitude_type), collapse = ", "), "\n")
sel <- sel[txt(receiver_type) == "eh"]
cat("  channel EH:", nrow(sel), "\n")
sel <- sel[txt(source_magnitude_type) == "md"]
cat("  magnitude MD:", nrow(sel), "\n")
sel <- sel[!is.na(num(p_arrival_sample))]
cat("  with P arrival:", nrow(sel), "\n")
if (nrow(sel) == 0) stop("no trace selected: see the counts above")
rm(meta); gc()
sel[, chunk := "seisbench"]
sel[, bucket := sub("\\$.*$", "", trace_name)]
sel[, row_in_bucket := as.integer(sub("^[^$]*\\$([0-9]+),.*$", "\\1", trace_name)) + 1L]
cat("selected traces:", nrow(sel), " events:", uniqueN(sel$source_id),
    " stations:", uniqueN(sel$receiver_code), "\n")

## ---- 2. waveform format ----------------------------------------------------------------
fid <- H5Fopen(h5_file, flags = "H5F_ACC_RDONLY")
read_format <- function(name, default) {
  out <- tryCatch(h5read(fid, paste0("data_format/", name)),
                  error = function(e) default)
  return(as.vector(out))
}
component_order <- as.character(read_format("component_order", "ZNE"))
sampling_rate   <- as.numeric(read_format("sampling_rate", fs))
cat("component order:", component_order, "  sampling rate:", sampling_rate, "\n")
if (abs(sampling_rate - fs) > 1e-6) stop("unexpected sampling rate: ", sampling_rate)
z_index <- regexpr("Z", component_order)[1]
if (z_index < 1) stop("no Z component in component_order = ", component_order)

read_vertical <- function(bucket, row) {
  # rhdf5 reverses the dimensions with respect to Python: a bucket of shape
  # (traces, 3, samples) is seen in R as (samples, 3, traces).
  arr <- tryCatch(h5read(fid, paste0("data/", bucket), index = list(NULL, NULL, row)),
                  error = function(e) NULL)
  if (is.null(arr)) {
    arr <- h5read(fid, paste0("data/", bucket), index = list(row, NULL, NULL))
  }
  w <- drop(arr)                               # samples x 3 or 3 x samples
  if (nrow(w) == 3) w <- t(w)
  return(as.numeric(w[, z_index]))
}

## ---- 3. preprocessing (as in app_00_extract_chunk.R) --------------------------------------
n_pre  <- round(pre_p_sec * fs)
n_post <- round(post_p_sec * fs)
n_win  <- n_pre + n_post
cosine_taper <- function(n, frac) {
  w <- rep(1, n)
  m <- floor(frac * n)
  if (m > 0) {
    ramp <- 0.5 * (1 - cos(pi * (0:(m - 1)) / m))
    w[1:m] <- ramp
    w[(n - m + 1):n] <- rev(ramp)
  }
  return(w)
}
taper <- cosine_taper(n_win, taper_frac)
if (have_signal) bp <- signal::butter(4, band_hz / (fs / 2), type = "pass")
preprocess <- function(z) {
  tt <- seq_along(z)
  z <- residuals(lm(z ~ tt))
  z <- z * taper
  if (have_signal) z <- signal::filtfilt(bp, z)
  return(as.numeric(z))
}

## ---- 4. read the selected traces -----------------------------------------------------
Y <- matrix(NA_real_, nrow = nrow(sel), ncol = n_win)
ok <- rep(FALSE, nrow(sel))
for (i in seq_len(nrow(sel))) {
  z <- read_vertical(sel$bucket[i], sel$row_in_bucket[i])
  p <- as.integer(round(sel$p_arrival_sample[i]))
  first <- p - n_pre + 1
  last  <- p + n_post
  if (first >= 1 && last <= length(z)) {
    Y[i, ] <- preprocess(z[first:last])
    ok[i] <- TRUE
  }
  if (i %% 1000 == 0) cat("  ", i, "of", nrow(sel), "\n")
}
H5Fclose(fid)

extract <- list(meta = sel[ok], Y = Y[ok, , drop = FALSE], chunk = "seisbench",
                fs = fs, pre_p_sec = pre_p_sec, post_p_sec = post_p_sec,
                band_hz = band_hz, filtered = have_signal)
saveRDS(extract, "results/stead_extract_seisbench.rds")
cat("saved results/stead_extract_seisbench.rds:", sum(ok), "traces\n")
