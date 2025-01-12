#!/bin/bash

chromosomes=('3' '17' '6' '10' '2' '12')
input_file="../inputs/gencode.vM25.annotation.gtf"

for chr in "${chromosomes[@]}"; do
    output_file="../inputs/gencode_chr${chr}_M25.gtf"
    grep "chr${chr}" "$input_file" > "$output_file"
    echo "Generated: $output_file"
done
