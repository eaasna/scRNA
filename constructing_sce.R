library(SingleCellExperiment)
library(SummarizedExperiment)

# decidua reference dataset
load(file="/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/decidua_5000_variably_expressed_counts.RData")
decidua_meta = read.table("/icgc/dkfzlsdf/analysis/B210/Evelin/E-MTAB-6701_arrayexpress_10x_meta.txt", header = TRUE)
decidua_meta = decidua_meta[which(decidua_meta$Cell %in% colnames(decidua_counts)), ]

library(tidyverse)
feature_symbol = str_split_fixed(row.names(decidua_counts), "_", 2)[,1]
rowData = data.frame(feature_symbol = feature_symbol, gene_id = str_split_fixed(row.names(decidua_counts), "_", 2)[,2])
row.names(decidua_counts) = feature_symbol

se = SummarizedExperiment(assays = list(logcounts = matrix(as.numeric(unlist(decidua_counts)),nrow=nrow(decidua_counts))), rowData = rowData, colData = decidua_meta)
rm(decidua_counts)

decidua = as(se, "SingleCellExperiment")
rm(se)
save(decidua, file = "/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/decidua_sce.RData")

# menstrual fluid dataset
library(Seurat)
load(file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/joined.RData"))

rowData = as.data.frame(row.names(seu[["RNA"]]))
colnames(rowData) <- c("feature_symbol")
colData = as.data.frame(colnames(seu[["RNA"]]))
colnames(colData) = c("Barcode")

seu_matrix = as.matrix(seu[["RNA"]]@counts)
se = SummarizedExperiment(assays = list(logcounts = matrix(as.numeric(unlist(seu_matrix)),nrow=nrow(seu_matrix))), rowData = rowData, colData = colData)
rm(seu)

menstrual = as(se, "SingleCellExperiment")
rm(se)
save(menstrual, file = "/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/menstrual_sce.RData")
