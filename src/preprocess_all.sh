#!/usr/bin/bash

set -e

# Define an array of file names
#files=("HRR1795848" "HRR1795849" "HRR1795850" "HRR1795851" "HRR1795852" "HRR1795853" "HRR1795854" "HRR1795855" "HRR1795856")
samples=("HRR1795848")

# Iterate over the array
for sample in "${samples[@]}"; do
    echo "Processing $sample"

    src/scLinaX_pipeline.sh /data/YH/Graves_dataset/CellRanger_results/${sample} ${sample}

    mkdir to_scp/${sample}

    samtools view -b -h /data/YH/Graves_dataset/CellRanger_results/${sample}/possorted_genome_X.sorted.bam chrX:79171491-79171491 -o to_scp/${sample}/${sample}_s162p_subset.bam

    cp /data/YH/Graves_dataset/CellRanger_results/${sample}/possorted_genome_X.sorted.bam to_scp/${sample}/${sample}_possorted_genome_X.sorted.bam

done
