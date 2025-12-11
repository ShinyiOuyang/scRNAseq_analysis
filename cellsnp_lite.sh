#!/usr/bin/bash

#scRNA-seq data info
barcode_rna="../HRR1795888/filtered_feature_bc_matrix/barcodes.tsv"
bam_rna="../HRR1795888/possorted_genome_X.sorted.bam"
UMITAG_rna="UB"

#run cellsnp for RNA
cellsnp-lite -s ${bam_rna} \
--minCOUNT 20 -b ${barcode_rna} \
--UMItag ${UMITAG_rna} -p 4 \
--minMAF 0.05 \
--minMAPQ 20 \
--refseq ../refdata-gex-GRCh38-2024-A/fasta/genome.fa \
--chrom=chrX -O RNA.chrX.snp.call

#process output for RNA
Rscript RCODE_process_cellsnp.r \
RNA.chrX.snp.call/cellSNP.samples.tsv \
RNA.chrX.snp.call/cellSNP.tag.AD.mtx \
RNA.chrX.snp.call/cellSNP.tag.DP.mtx \
RNA.chrX.snp.call/cellSNP.tag.OTH.mtx \
RNA.chrX.snp.call/cellSNP.base.vcf \
RNA.chrX.snp.call.summary.tsv.gz