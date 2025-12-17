#!/usr/bin/bash

set -e

# Following steps from https://ytomofuji.github.io/scLinaX/articles/scLinaX_preprocessing_example.html

# Step 1: Data Preparation
#python3 ./src/get_chrX_reads.py $1 # Extract chrX reads from BAM

# Unzip filtered feature barcode matrix
if [ ! -d "$1/filtered_feature_bc_matrix/" ]; then 
    mkdir $1/filtered_feature_bc_matrix/
fi
if [ ! -e "$1/filtered_feature_bc_matrix.tar.gz" ]; then
    tar -xf $1/filtered_feature_bc_matrix.tar.gz -C $1/filtered_feature_bc_matrix/
fi
if [ ! -e "$1/filtered_feature_bc_matrix/barcodes.tsv.gz" ]; then
    zcat $1/filtered_feature_bc_matrix/barcodes.tsv.gz > $1/filtered_feature_bc_matrix/barcodes.tsv
fi

# Step 2: Run Cellsnp-lite
bash ./src/cellsnp_lite.sh $1 $2 

# Step 3: Run Annovar
bash ./src/annovar.sh

# Step 4: QC
Rscript ./src/RCODE_qc.R RNA.chrX.snp.call.summary.tsv.gz RNA.chrX.snp.call.hg38_multianno.txt $3