library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
#library(ArchR)
library(harmony)
library(ggpubr)

graves1.data <- Read10X(data.dir = "/20251029_Graves/filtered_feature_bc_matrix/")
graves1 <- CreateSeuratObject(counts = E_H1, project = "E_H-01", min.cells = 3, min.features = 400)

