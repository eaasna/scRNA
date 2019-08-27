#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


norm = "SCT"
#norm = args[1]


library(Seurat)
library(devtools)

path = "/icgc/dkfzlsdf/analysis/B210/Evelin/decidua/"
seu <- Read10X(path)
seu <- CreateSeuratObject( seu )

seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")
high_mt_cells = names(seu$nFeature_RNA[seu$percent.mt > 20])
seu = subset(seu, cells = high_mt_cells, invert = TRUE)

# seurat slot for cell type annotation
metadata = read.table("/icgc/dkfzlsdf/analysis/B210/Evelin/E-MTAB-6701_arrayexpress_10x_meta.txt", header = TRUE)
metadata = metadata[which(metadata$Cell %in% colnames(seu[["RNA"]])), ]
seu$annotation = metadata$annotation


if (norm == "log"){
  seu <- NormalizeData(seu, normalization.method = "LogNormalize" )
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, assay = "RNA")
  seu <- ScaleData(seu, features = rownames(seu) )
  nPCA = 20
  
}
if (norm == "SCT"){
  # size limit: 2048 * 1024^2 
  options(future.globals.maxSize=891289600)
  seu <- SCTransform(seu, verbose = FALSE, conserve.memory = TRUE )
  nPCA = 30
}


seu <- RunPCA(seu, features = VariableFeatures(seu) )
#seu <- RunUMAP(seu, dims = 1:nPCA)
seu <- RunTSNE(seu, dims = 1:nPCA )

seu <- FindNeighbors(seu, dims = 1:nPCA)
seu <- FindClusters(seu, verbose = FALSE)
save( seu, file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/decidua/",norm,"_seu.RData" ))