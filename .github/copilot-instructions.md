# scRNAseq X Chromosome Inactivation (XCI) Analysis - AI Agent Instructions

## Project Overview
This is a bioinformatics pipeline analyzing X chromosome inactivation (XCI) escape using single-cell RNA-seq data. The project integrates three components:
1. **scLinaX R package** (`scLinaX/`) - Core XCI escape quantification library
2. **Analysis pipeline** (`src/`) - Multi-stage processing: quality control → allele-specific expression → XCI analysis
3. **Analysis scripts** (Rmd files, R scripts in root) - Sample-specific analyses and visualization

## Key Architecture Patterns

### Data Flow
```
Raw BAM/HDF5 reads 
  → cellsnp-lite (SNP genotyping per cell) 
  → ANNOVAR (variant annotation) 
  → QC filtering (RCODE_qc.R) 
  → scLinaX analysis (run_scLinaX) 
  → Seurat/Azimuth (cell type annotation) 
  → DEG analysis
```

### Directory Structure
- `input/` - Sample folders with CellRanger outputs (HDF5 matrices, TSV SNP calls)
- `output/` - Sample-specific results (CSVs from scLinaX, plots, cell type annotations)
- `src/` - Bash pipelines and utility R scripts (not meant for direct execution; see `scLinaX_pipeline.sh`)
- `scLinaX/` - R package with core functions (`run_scLinaX`, `summarize_scLinaX`, QC utilities)

### Sample Processing Convention
Each sample (e.g., HRR1795848) has:
- Input: `input/SAMPLE/filtered_feature_bc_matrix.h5` + `input/SAMPLE/SAMPLE_RNA_QC_passed_SNP_df.tsv`
- Outputs: `output/SAMPLE/` with CSVs (plain results, cell counts, SNP summary, phasing, clustering)
- Pattern: Loop through sample list in analysis scripts (e.g., `src/scLinaX_analysis.R`)

## Critical Dependencies & Workflows

### R Package Dependencies (via renv)
- **scLinaX** (custom, local) - XCI escape quantification; implements `run_scLinaX()`, `summarize_scLinaX()`, QC functions
- **Seurat/Azimuth** - scRNA-seq analysis, cell type annotation (uses pre-trained PBMC reference)
- **Tidyverse stack** (dplyr, readr, stringr, purrr, magrittr) - Data wrangling
- **igraph** - Network clustering for allele phasing

### Python/Bash Tools
- `cellsnp-lite` - Multi-sample SNP calling per cell (installed via conda)
- `ANNOVAR` - Variant annotation (must be in `PATH`)
- `samtools` - BAM manipulation (reads, indexing)
- `FastQC` - Quality control (optional)

### Environment Setup (One-time)
```bash
# Set conda environment
source activate cellsnp_lite

# In R: restore renv (syncs to lockfile)
renv::restore()
```

## scLinaX Package API (Core Functions)

### Main Analysis
- `run_scLinaX(ASE_df, XCI_ref, QCREF, ...)` - Core XCI analysis
  - Requires: allele-specific expression dataframe with columns: `SNP_ID`, `Sample_ID`, `REFcount`, `ALTcount`, `cell_barcode`, `POS`, `Gene`
  - Key params: `SNP_DETECTION_DP=30`, `SNP_DETECTION_MAF=0.1`, `HE_allele_cell_number_THR=50`, `PVAL_THR=0.01`, `RHO_THR=0.5`
  - Returns: List with `$result` (per-gene-cell assignments), `$Max_Num_Table_result` (cell counts), `$phasing_result` (allele phasing), `$clustering_result` (K-means clusters)

- `summarize_scLinaX(scLinaX_res, QC_total_allele_THR=0, Annotation=NULL)` - Aggregate results by cell type
  - Input: scLinaX result object + optional cell type annotation dataframe
  - Output: Dataframe with gene × cell type summary including `Gene_class` (PAR1/nonPAR_escape/nonPAR_inactive/etc.)

### Reference Data (Built-in)
- `data("XCI_ref")` - Gene-level XCI status (escape/inactive/variable/unknown)
- `data("AIDA_QCREF")` - QC reference for removing potential escapee genes in reference set
- Modify XCI_ref status if needed (e.g., `XCI_ref$XCI_status[XCI_ref$Gene == "GENE_NAME"] <- "escape"`)

## Project-Specific Patterns

### QC Strategy
1. **cellsnp-lite QC** - Filters by read depth (DP≥30), minor allele frequency (MAF 0.1-0.9)
2. **ANNOVAR QC** - Retains intronic/UTR/exonic variants; removes multi-gene SNPs
3. **scLinaX internal QC** - `make_QCed_df()` applies allele ratio thresholds; `Gene_QC()` removes genes with low Xi expression
4. **scLinaX output QC** - Filter by `QC_total_allele_THR` when summarizing

### Visualization Conventions
- Gene classification boxplots: X-axis ordered as `PAR1 → nonPAR_escape → nonPAR_variable → nonPAR_inactive → nonPAR_unknown → PAR2`
- Color scheme (from `scLinaX_analysis.R`): PAR1=#ae2d68, escape=#ff1700, variable=#ffa600, inactive=#9db300, unknown=#666666, PAR2=#729efd
- Allele barchart: Xa (paternal/maternal assignment) grouped by frequency
- Save plots: `ggsave(filename, device="png", width=6, height=4)`

### Cell Type Annotation Workflow
If Azimuth prediction not pre-computed:
```r
Prediction <- RunAzimuth(seurat_obj, "pbmcref", assay = "RNA")
CellTypes_simple <- Prediction@meta.data[, c("cell_barcode", "predicted.celltype.l2")]
# Save as SAMPLE_CellTypes_Azimuth_simple.csv
```

## Common Tasks & Patterns

### Running Full Pipeline for a Sample
```bash
bash src/scLinaX_pipeline.sh /data/path/to/sample SAMPLE_NAME
```
Automatically: extracts chrX reads → runs cellsnp → ANNOVAR → QC → outputs `SAMPLE_RNA_QC_passed_SNP_df.tsv`

### Iterating Over Samples
Define sample vector in R: `samples <- c('HRR1795848', 'HRR1795849', ...)`
Loop through with dynamic path construction: `curr_input_dir = paste0(input_dir, samplename, "/")`

### Merging Seurat with X Chromosome Annotations
```r
data <- Read10X_h5(paste0(input_dir, samplename, "/filtered_feature_bc_matrix.h5"))
data <- CreateSeuratObject(data, project = "ProjectName")
xchrom_annotation <- read.csv(paste0(output_dir, samplename, ".xchrom_annotation.csv"))
data <- merge(data, xchrom_annotation, by = "cell_barcode")  # Adds Xa column to metadata
```

### Standardizing Output CSVs
Always include: `cell_barcode`, relevant allele assignments, gene metadata
Use `write.csv(..., row.names = FALSE)` to ensure readability
Save to `output/SAMPLE/` with naming pattern: `SAMPLE.ANALYSIS_TYPE.csv` (e.g., `HRR1795848.sclinax_plain_results.csv`)

## Development Notes
- **Custom package**: scLinaX source in `scLinaX/R/` - modify functions directly, then reinstall with `devtools::load_all()` or rebuild
- **renv cache**: Check with `renv::paths$cache()` if packages not installing
- **Memory management**: Large HDF5 matrices loaded into memory; monitor with `object.size()` for >10GB samples
- **Reproducibility**: All analysis scripts should declare sample list explicitly; avoid hardcoded indices
