#!/usr/bin/env bash

for var in *.bw; do
    if [ -e "$(basename "$var" .bw).bedgraph" ]; then
        continue
    fi
    bigWigToBedGraph "$var" "$(basename "$var" .bw).bedgraph"
done
