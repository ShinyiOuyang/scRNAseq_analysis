#!/usr/bin/bash

set -e

# Define an array of file names
files=("HRR1795848" "HRR1795849" "HRR1795850" "HRR1795851" "HRR1795852" "HRR1795853", "HRR1795854", "HRR1795855", "HRR1795856")

# Iterate over the array
for file in "${files[@]}"; do
    echo "Processing $file"

    echo src/scLinaX_pipeline.sh /data/YH/Graves_dataset/CellRanger_results_SO_copy/${file} ${file}
done
