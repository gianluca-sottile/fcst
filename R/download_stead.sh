#!/bin/bash
# =============================================================================
#  download_stead.sh
#  Downloads one STEAD chunk at a time, extracts the traces needed for the
#  paper (app_00_extract_chunk.R) and then deletes the chunk to save space.
#
#  Usage (from the R/ folder, in the Terminal):
#     bash download_stead.sh            # all chunks 2-6
#     bash download_stead.sh chunk3     # a single chunk
#
#  Space needed: about 30 GB free at a time (one zipped chunk + its content).
#  Time: each chunk is about 14 GB to download.
# =============================================================================
set -e
DATA_DIR="$HOME/Data/STEAD"          # must match data_dir in app_00_extract_chunk.R
KEEP_CHUNKS="no"                     # "yes" to keep the original files
CHUNKS=${1:-"chunk2 chunk3 chunk4 chunk5 chunk6"}

mkdir -p "$DATA_DIR"
for CH in $CHUNKS; do
  echo "=== $CH ==="
  if [ ! -f "results/stead_extract_$CH.rds" ]; then
    ZIP="$DATA_DIR/$CH.zip"
    if [ ! -f "$DATA_DIR/$CH/$CH.hdf5" ]; then
      # -L follows the redirect of the short link; -C - resumes interrupted downloads
      curl -L -C - -o "$ZIP" "https://rebrand.ly/$CH"
      SIZE=$(stat -f%z "$ZIP" 2>/dev/null || stat -c%s "$ZIP")
      if [ "$SIZE" -lt 1000000000 ]; then
        echo "The downloaded file is too small ($SIZE bytes): the link probably"
        echo "opened a web page. Download $CH from the browser (https://rebrand.ly/$CH),"
        echo "put it in $DATA_DIR and run this script again."
        exit 1
      fi
      unzip -o "$ZIP" -d "$DATA_DIR"
      rm -f "$ZIP"
    fi
    # the archive may unpack to $DATA_DIR/$CH/ or directly to $DATA_DIR/
    if [ -f "$DATA_DIR/$CH.hdf5" ]; then
      mkdir -p "$DATA_DIR/$CH"
      mv "$DATA_DIR/$CH.hdf5" "$DATA_DIR/$CH.csv" "$DATA_DIR/$CH/"
    fi
    Rscript app_00_extract_chunk.R "$CH"
  else
    echo "extract already present, skipped"
  fi
  if [ "$KEEP_CHUNKS" != "yes" ]; then
    rm -rf "$DATA_DIR/$CH"
  fi
done
echo "Done. Now run app_01_prepare_stead.R and app_02_fcst_stead.R."
