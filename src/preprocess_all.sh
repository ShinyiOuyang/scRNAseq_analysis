#!/usr/bin/bash

set -e

bash src/scLinaX_pipeline.sh /data/CellRanger_results_SO_copy/HRR1795892 /data/YH/Graves_dataset/CellRanger_results_SO_copy/refdata-gex-GRCh38-2024-A/fasta/genome.fa HRR1795892

bash src/scLinaX_pipeline.sh /data/CellRanger_results_SO_copy/HRR1795893 /data/YH/Graves_dataset/CellRanger_results_SO_copy/refdata-gex-GRCh38-2024-A/fasta/genome.fa HRR1795893

bash src/scLinaX_pipeline.sh /data/CellRanger_results_SO_copy/HRR1795894 /data/YH/Graves_dataset/CellRanger_results_SO_copy/refdata-gex-GRCh38-2024-A/fasta/genome.fa HRR1795894

bash src/scLinaX_pipeline.sh /data/CellRanger_results_SO_copy/HRR1795895 /data/YH/Graves_dataset/CellRanger_results_SO_copy/refdata-gex-GRCh38-2024-A/fasta/genome.fa HRR1795895