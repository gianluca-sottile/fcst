################################################################################
#  app_00_extract_chunk.R
#
#  Step 0 of the STEAD application: extracts, from ONE downloaded chunk, the
#  traces that satisfy the selection of the paper, and saves only what is
#  needed (metadata + preprocessed vertical window) in a small RDS file.
#  After this script has run on a chunk, the chunk (~14 GB) can be deleted.
#
#  Usage: set 'chunk' below (or pass it on the command line:
#         Rscript app_00_extract_chunk.R chunk2) and run once per chunk,
#         chunk2, ..., chunk6 (chunk1 is noise only).
#
#  Input : <data_dir>/<chunk>/<chunk>.csv and <data_dir>/<chunk>/<chunk>.hdf5
#  Output: results/stead_extract_<chunk>.rds
#
#  Requirements: data.table, rhdf5 (BiocManager::install("rhdf5")),
#                signal (optional, band-pass filter)
################################################################################

library(data.table)
library(rhdf5)
have_signal <- requireNamespace("signal", quietly = TRUE)

## ---- settings ----------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
chunk    <- if (length(args) > 0) args[1] else "chunk2"
data_dir <- "~/Data/STEAD"

lat_range  <- c(45, 47)
lon_range  <- c(-123, -121)
fs         <- 100      # Hz
pre_p_sec  <- 0.5      # window: 0.5 s before P ...
post_p_sec <- 12       # ... to 12 s after P
band_hz    <- c(1, 20)
taper_frac <- 0.05

dir.create("results", showWarnings = FALSE)

## ---- metadata and selection ---------------------------------------------------------
csv_file <- path.expand(file.path(data_dir, chunk, paste0(chunk, ".csv")))
h5_file  <- path.expand(file.path(data_dir, chunk, paste0(chunk, ".hdf5")))
if (!file.exists(csv_file) || !file.exists(h5_file)) {
  stop("files not found for ", chunk, ": check data_dir and the folder names")
}
meta <- fread(csv_file)
sel <- meta[
  trace_category == "earthquake_local" &
  source_latitude  >= lat_range[1] & source_latitude  <= lat_range[2] &
  source_longitude >= lon_range[1] & source_longitude <= lon_range[2] &
  receiver_type == "EH" &
  source_magnitude_type == "MD" &
  !is.na(p_arrival_sample)
]
sel[, chunk := chunk]
cat(chunk, ": selected traces", nrow(sel), " events", uniqueN(sel$source_id), "\n")

## ---- preprocessing helpers ----------------------------------------------------------
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
  z <- residuals(lm(z ~ tt))            # mean and linear trend removed
  z <- z * taper
  if (have_signal) z <- signal::filtfilt(bp, z)
  return(as.numeric(z))
}

read_vertical <- function(fid, trace_name) {
  w <- h5read(fid, paste0("data/", trace_name))
  if (nrow(w) == 3) w <- t(w)           # rhdf5 may return 3 x 6000
  return(as.numeric(w[, 3]))            # columns: E, N, Z
}

## ---- waveforms ---------------------------------------------------------------------
Y <- matrix(NA_real_, nrow = nrow(sel), ncol = n_win)
ok <- rep(FALSE, nrow(sel))
if (nrow(sel) > 0) {
  fid <- H5Fopen(h5_file, flags = "H5F_ACC_RDONLY")
  for (i in 1:nrow(sel)) {
    z <- read_vertical(fid, sel$trace_name[i])
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
}
extract <- list(meta = sel[ok], Y = Y[ok, , drop = FALSE], chunk = chunk,
                fs = fs, pre_p_sec = pre_p_sec, post_p_sec = post_p_sec,
                band_hz = band_hz, filtered = have_signal)
out_file <- file.path("results", paste0("stead_extract_", chunk, ".rds"))
saveRDS(extract, out_file)
cat("saved", out_file, ":", sum(ok), "traces\n")
