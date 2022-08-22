#!/usr/bin/env bash

# conda activate anaconda

cell_types=("Zygote" "2cell" "4cell" "8cell" "Morula" "ICM" "TE")
for each_type in "${cell_types[@]}"; do
    echo "Processing ${each_type}"
    python NOMe_bulk_test.py "${each_type}_merged_gch.bedgraph" > "${each_type}_special_peaks.bedgraph"
    lines=$(awk 'END { print NR }' "${each_type}_special_peaks.bedgraph")
    echo "The output ${each_type} file has ${lines} lines"
done
