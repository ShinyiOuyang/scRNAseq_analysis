# library(remotes)
library(Seurat)
library(Azimuth)
library(SeuratData)

available_data <- AvailableData()
options(timeout = 1000)
#InstallData("pbmcref")

# setwd("/Users/yun-huang/Documents/GPR174/scRNAseq/20250712_resequence_analysis")
# samples = c("C2_resequence","C3_resequence","HD89_resequence")
#input_dir <- "/data/YH/Graves_dataset/CellRanger_results/"
samples <- c('HRR1795846')
input_dir <- "/Users/shinyiouyang/Downloads/scRNAseq_analysis/input/"
output_dir <- "/Users/shinyiouyang/Downloads/scRNAseq_analysis/output/"

for (samplename in samples) {
    curr_output_dir = paste0(output_dir, samplename, "/")
    data <- Read10X_h5(paste0(input_dir, samplename, "/filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
    data <- CreateSeuratObject(data, project = "Graves_scRNAseq", assay = "RNA")
    xchrom_annotation = read.csv(paste0(curr_output_dir, samplename, ".xchrom_annotation.csv"), header = TRUE)

    Prediction <- RunAzimuth(data, "pbmcref", assay = "RNA")

    CellTypes_Azimuth <- Prediction@meta.data
    CellTypes_Azimuth$Barcode <- rownames(CellTypes_Azimuth)
    CellTypes_Azimuth_simple <- CellTypes_Azimuth[, c(12, 8)]

    colnames(CellTypes_Azimuth_simple) <- c("cell_barcode", "predicted.celltype.l2")
    data$predicted.celltype.l2 <- CellTypes_Azimuth$predicted.celltype.l2[match(colnames(data), CellTypes_Azimuth$Barcode)]

    cell_barcodes <- colnames(data)

    barcode_mapping <- data.frame(
        cell_barcode = cell_barcodes,
        stringsAsFactors = FALSE
    )

    barcode_mapping <- barcode_mapping %>%
    left_join(xchrom_annotation, by = "cell_barcode") %>%
    mutate(
        # Replace NA with "Uncategorized" for cells without annotation
        Xa = replace_na(Xa, "Uncategorized")  # adjust column name to match yours
    )

    data$Xa <- barcode_mapping$Xa[match(colnames(data), barcode_mapping$cell_barcode)]

    #data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
    #qc_plot = VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    #plot(qc_plot)
    data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

    # Subset to NK cells
    seurat_nk <- subset(data, predicted.celltype.l2 == "NK")

    # Set Xa_status as the active identity
    Idents(seurat_nk) <- "Xa"

    # Run FindMarkers: Allele A vs Allele B
    nk_deg <- FindMarkers(
    seurat_nk,
    ident.1 = "Allele_A",
    ident.2 = "Allele_B",
    min.pct = 0.1,        # Express in ≥10% of cells in either group
    logfc.threshold = 0.25
    )

    # Add gene names and filter
    nk_deg$gene <- rownames(nk_deg)
    nk_deg_sig <- nk_deg %>% filter(p_val_adj < 0.05)
    print(head(nk_deg))
    print(head(nk_deg_sig))

    write.csv(nk_deg, paste0(curr_output_dir, samplename, "_NK_deg.csv"), row.names = FALSE)
}