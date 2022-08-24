#!/usr/bin/env bash
bedtools merge -i ICM_ATAC_peaks.bed > tmpfile && mv tmpfile ICM_ATAC_peaks.bed;
printf "ICM done\n"
bedtools merge -i TE_ATAC_peaks.bed > tmpfile && mv tmpfile TE_ATAC_peaks.bed;
printf "TE done\n"


