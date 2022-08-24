#!/usr/bin/env bash
touch ICM_ATAC_peaks.bed
cat *ICM*.bed | sort -u -k1,1 -k2n > ICM_ATAC_peaks.bed

touch TE_ATAC_peaks.bed
cat *TE*.bed | sort -u -k1,1 -k2n > TE_ATAC_peaks.bed

