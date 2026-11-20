## Installs the packages used by the scripts (the core of FCST uses base R only).
cran <- c("mclust", "cluster", "tclust", "clustEff", "data.table", "signal", "maps", "GET")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) install.packages(missing)
# rhdf5 (Bioconductor) is needed only to read the STEAD waveforms
if (!requireNamespace("rhdf5", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("rhdf5")
}
