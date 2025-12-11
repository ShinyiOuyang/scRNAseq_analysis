#!/usr/bin/bash

annovarpath=/home/hannah/annovar/

#  scRNA-seq
zcat RNA.chrX.snp.call.summary.tsv.gz | \
awk -F"\t" 'NR>1 {print $2"\t"$3"\t"$3"\t"$4"\t"$5"\t"$1}' | \
LANG=C sort | LANG=C uniq \
> RNA.chrX.snp.call.summary.annov.input

$annovarpath/table_annovar.pl \
RNA.chrX.snp.call.summary.annov.input \
$annovarpath/humandb/ -buildver hg38 \
-dbtype refGeneWithVer \
-out RNA.chrX.snp.call \
-remove -protocol refGene,ensGene,avsnp150,ALL.sites.2015_08 \
-operation g,g,f,f -nastring NA -polish

# # scATAC-seq
# zcat ATAC.chrX.snp.call.summary.tsv.gz | \
# awk -F"\t" 'NR>1 {print $2"\t"$3"\t"$3"\t"$4"\t"$5"\t"$1}' | \
# LANG=C sort | LANG=C uniq \
# > ATAC.chrX.snp.call.summary.annov.input
# /path/to/annovar/table_annovar.pl \
# ATAC.chrX.snp.call.summary.annov.input \
# /path/to/annovar/humandb/ -buildver hg38 \
# -out ATAC.chrX.snp.call \
# -remove -protocol refGene,ensGene,avsnp150,ALL.sites.2015_08 \
# -operation g,g,f,f -nastring NA -polish