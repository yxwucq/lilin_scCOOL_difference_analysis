#!/usr/bin/env bash

# 将每个类型输出的bedgraph file合并在一起
cell_types=("Zygote" "2cell" "4cell" "8cell" "Morula" "ICM" "TE")

for each_type in "${cell_types[@]}"; do

    #create file
    echo "Processing $each_type"
    file_name="${each_type}_merged_gch.bedgraph"
    if [ -e "$file_name" ]; then #check&remove
        rm "$file_name"
    fi
    touch "$file_name"

    # cat together
    for each_type_file in *"_$each_type"*.bedgraph; do
        awk '$1 ~ /chr[0-9]/ { print }' "$each_type_file" >> "$file_name" #去除chrX,chrY,lambda,MIT
    done
    echo "Processing $each_type: Merged"
    # sort&new_columns
    sort -k1,1 -k2n "$file_name" > tmpfile && mv tmpfile "$file_name" #sort the key
    echo "Processing $each_type: Sorted"

    # count the reads of samples
    awk '{print $1, $2, $3, 1-$4, $4} 
        { count_un_methyl += (1-$4); count_methyl += $4}
        END { print count_un_methyl ; print count_methyl }
        ' "$file_name" > tmpfile && mv tmpfile "$file_name"
    echo "Processing $each_type: Total Reads Counts Finished"
    count_array=( $( tail -2 < "$file_name" ) )
    echo "${count_array[@]}" > "${each_type}_merged_gch_total_count.xls" # save the count number to attached file
    head -n -2 "$file_name" > tmpfile && mv tmpfile "$file_name" # cut the last two lines off

done