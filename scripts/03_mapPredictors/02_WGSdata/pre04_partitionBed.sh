#!/bin/bash

# Set base directory where your input files are located
base_dir="data/MutTables/WholeGenomeData/"
cd "$base_dir" || exit

# Loop through chromosomes 1 to 22
for chr in {1..22}; do
    infile="MutsResult_chr${chr}.bed"

    # Check if the file exists
    if [[ ! -f "$infile" ]]; then
        echo "File $infile not found!"
        continue
    fi

    # Split into exactly 100 parts
    split -d -n l/100 --additional-suffix=".bed" "$infile" "genomeMuts_chr${chr}_part"

    # Rename files with proper numbering
    # Use sort -V for natural sort of numbers in strings
    part_num=1
    for f in $(ls genomeMuts_chr${chr}_part*.bed | sort -V); do
        new_name=$(printf "genomeMuts_chr%s_part%03d.bed" "$chr" "$part_num")
        mv "$f" "$new_name"
        ((part_num++))
    done
done

# Move output files to temp2 directory
mkdir -p temp2
mv genomeMuts* ./temp2