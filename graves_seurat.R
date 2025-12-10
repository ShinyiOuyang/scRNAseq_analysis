library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
#library(ArchR)
library(harmony)
library(ggpubr)

graves1.data <- Read10X(data.dir = "../HRR1795847/filtered_feature_bc_matrix/")
graves1 <- CreateSeuratObject(counts = graves1.data, project = "Graves_1", min.cells = 3, min.features = 400)
graves1[["percent.mt"]] <- PercentageFeatureSet(graves1, pattern = "^MT-") 
graves1 <- subset(graves1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)

DefaultAssay(graves1) <- 'RNA'

graves1 <- NormalizeData(graves1, normalization.method = "LogNormalize", scale.factor= 1e4) # default parameters for NormalizeData
graves1 <- FindVariableFeatures(graves1, nfeatures = 4000) # default is nfeatures = 2000                     
graves1 <- ScaleData(graves1, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
graves1 <- RunPCA(graves1, npcs = 50, verbose = FALSE)
graves1 <- scRNA %>% RunHarmony("orig.ident", plot_convergence = TRUE)
graves1 <- RunUMAP(graves1, reduction = "harmony", dims = 1:30)
graves1 <- FindNeighbors(graves1, reduction = "harmony", dims = 1:50)
graves1 <- FindClusters(graves1, resolution = 0.6)