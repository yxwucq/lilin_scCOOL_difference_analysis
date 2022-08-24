#! /bin/sh

touch ICM_special_peaks.bed
touch TE_special_peaks.bed

cp ICM_ATAC_peaks.bed ICM_special_peaks.bed
cp TE_ATAC_peaks.bed TE_special_peaks.bed

tail ICM_special_peaks.bed

for variable in ./*Zygote*.bed;
do bedtools intersect -v -a ICM_special_peaks.bed -b "$variable" > tmpfile && mv tmpfile ICM_special_peaks.bed;
bedtools intersect -v -a TE_special_peaks.bed -b "$variable" > tmpfile && mv tmpfile TE_special_peaks.bed;
done

echo "Zygote done"

for variable in ./*2cell*.bed;
do bedtools intersect -v -a ICM_special_peaks.bed -b "$variable" > tmpfile && mv tmpfile ICM_special_peaks.bed;
bedtools intersect -v -a TE_special_peaks.bed -b "$variable" > tmpfile && mv tmpfile TE_special_peaks.bed;
done

echo "2cell done"

for variable in ./*4cell*.bed;
do bedtools intersect -v -a ICM_special_peaks.bed -b "$variable" > tmpfile && mv tmpfile ICM_special_peaks.bed;
bedtools intersect -v -a TE_special_peaks.bed -b "$variable" > tmpfile && mv tmpfile TE_special_peaks.bed;
done

echo "4cell done"

for variable in ./*8cell*.bed;
do bedtools intersect -v -a ICM_special_peaks.bed -b "$variable" > tmpfile && mv tmpfile ICM_special_peaks.bed;
bedtools intersect -v -a TE_special_peaks.bed -b "$variable" > tmpfile && mv tmpfile TE_special_peaks.bed;
done

echo "8cell done"

for variable in ./*Morula*.bed;
do bedtools intersect -v -a ICM_special_peaks.bed -b "$variable" > tmpfile && mv tmpfile ICM_special_peaks.bed;
bedtools intersect -v -a TE_special_peaks.bed -b "$variable" > tmpfile && mv tmpfile TE_special_peaks.bed;
done

echo "Morula done"
echo "All done"
