#!/usr/bin/env bash
# Do Not Directly Run This Code
# 吴宇轩 2022.8.19

# 这是对处理过程的代码总结，供复现和提供思路作用
# Title: Single-cell multi-omics sequencing of human early embryos
# Data Type: Methylation profiling by high throughput sequencing
# Data Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100272

# Pre-processing
## Unzip the files
# 选取其中代表可及性的GCT数据
tar -xf GSE100272_RAW.tar.1 --wildcards "*.GCT.bw" --exclude="*alpha*"

## converting the files to .bedgraph
# 利用bigWigToBedGraph转换文件格式，需要创建clean环境
# 记得切换环境
for var in *.bw; do
    if [ ! -e "$(basename "$var" .bw).bedgraph" ]; then
        continue
    fi
    bigWigToBedGraph "$var" "$(basename "$var" .bw).bedgraph"
done
# Monitoring the progress
while true
do 
    numbw=$(ls *.bw | wc -l)
    numbg=$(ls *.bedgraph | wc -l)
    if (($numbw==$numbg)); then
        break
    fi
    printf "Processing %d/%d\n" "$numbg" "$numbw"
    sleep 30
done

# 将每个类型输出的bedgraph file合并在一起
# 合并、消除X，Y，Mit、排序、计算GC甲基化比例
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
    for each_type_file in *"$each_type"*.bedgraph; do
        awk '$1 ~ /chr[0-9]/ { print }' "$each_type_file" >> "$file_name" #去除chrX,chrY,lambda,MIT
    done
    # sort&new_columns
    sort -k1,1 -k2n "$file_name" > tmpfile && mv tmpfile "$file_name" #sort the key
    echo "Processing $each_type: Merged"

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

# monitoring the process
while true; do
    sleep 60
    tail merge_for_peaks_calling.sh.o3339068 
done

