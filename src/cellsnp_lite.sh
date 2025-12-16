#!/usr/bin/bash

#scRNA-seq data info
barcode_rna="${1}/filtered_feature_bc_matrix/barcodes.tsv"
bam_rna="${1}/possorted_genome_X.sorted.bam"
UMITAG_rna="UB"

#run cellsnp in mode 2b
cellsnp-lite -s ${bam_rna} \
--minCOUNT 20 \
--cellTAG None \
--UMItag None \
-p 16 \
--minMAF 0.05 \
--minMAPQ 20 \
--refseq ${2} \
--chrom=chrX -O bulk.chrX.snp.call

#../refdata-gex-GRCh38-2024-A/fasta/genome.fa

#run cellsnp in mode 1a
cellsnp-lite -s ${bam_rna} \
--minCOUNT 20 -b ${barcode_rna} \
-R bulk.chrX.snp.call/cellSNP.base.vcf \
--UMItag ${UMITAG_rna} \
-p 16 \
--minMAF 0.05 \
--minMAPQ 20 \
--refseq ${2} \
--chrom=chrX -O RNA.chrX.snp.call


#process output for RNA
Rscript ./src/RCODE_process_cellsnp.r \
RNA.chrX.snp.call/cellSNP.samples.tsv \
RNA.chrX.snp.call/cellSNP.tag.AD.mtx \
RNA.chrX.snp.call/cellSNP.tag.DP.mtx \
RNA.chrX.snp.call/cellSNP.tag.OTH.mtx \
RNA.chrX.snp.call/cellSNP.base.vcf \
RNA.chrX.snp.call.summary.tsv.gz