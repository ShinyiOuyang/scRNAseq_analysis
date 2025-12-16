#!/usr/bin/bash

set -e

# Following steps from https://ytomofuji.github.io/scLinaX/articles/scLinaX_preprocessing_example.html

# Step 1: Data Preparation
#python3 ./get_chrX_reads.py $1 # Extract chrX reads from BAM

# Unzip filtered feature barcode matrix
mkdir $1/filtered_feature_bc_matrix/
tar -xf $1/filtered_feature_bc_matrix.tar.gz -C $1/filtered_feature_bc_matrix/
zcat $1/filtered_feature_bc_matrix/barcodes.tsv.gz > $1/filtered_feature_bc_matrix/barcodes.tsv

# Step 2: Run Cellsnp-lite
bash ./cellsnp_lite.sh $1

# Step 3: Run Annovar
bash ./annovar.sh

# Step 4: QC
Rscript RCODE_qc.R $1

