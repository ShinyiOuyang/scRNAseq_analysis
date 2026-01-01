# library(remotes)
library(Seurat)
library(Azimuth)
library(SeuratData)

available_data <- AvailableData()
options(timeout = 1000)
InstallData("pbmcref")

# setwd("/Users/yun-huang/Documents/GPR174/scRNAseq/20250712_resequence_analysis")
# samples = c("C2_resequence","C3_resequence","HD89_resequence")
#input_dir <- "/data/YH/Graves_dataset/CellRanger_results/"
input_dir <- "/Users/shinyiouyang/Downloads/scRNAseq_analysis/input/"
output_dir <- "/Users/shinyiouyang/Downloads/scRNAseq_analysis/output/"
samples <- c('HRR1795848','HRR1795849', 'HRR1795850', 'HRR1795855', 'HRR1795856')

for (samplename in samples) {
    curr_output_dir = paste0(output_dir, samplename, "/")
    data <- Read10X_h5(paste0(input_dir, samplename, "/filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
    data <- CreateSeuratObject(data, project = "Graves_scRNAseq", assay = "RNA")
    # data = SCTransform(data, assay = "RNA", verbose = FALSE)

    Prediction <- RunAzimuth(data, "pbmcref", assay = "RNA")

    CellTypes_Azimuth <- Prediction@meta.data
    CellTypes_Azimuth$Barcode <- rownames(CellTypes_Azimuth)
    CellTypes_Azimuth_simple <- CellTypes_Azimuth[, c(12, 8)]

    colnames(CellTypes_Azimuth_simple) <- c("cell_barcode", "Annotation")
    rownames(CellTypes_Azimuth_simple) <- NULL
    write.csv(CellTypes_Azimuth_simple, paste0(curr_output_dir, samplename, "_CellTypes_Azimuth_simple.csv"), row.names = FALSE)

    p1 <- DimPlot(Prediction, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend() + ggplot2::labs(title = paste0(samplename, " Azimuth Cell Type L2"))
    plot(p1)
}