#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(SingleCellExperiment)
library(SummarizedExperiment)
library(Seurat)
library(dplyr)

path = "/icgc/dkfzlsdf/analysis/B210/Evelin/"

origin = args[1]
type = args[2]
assay = args[3]

if (origin == "menstrual"){
  load(file = paste0(path, "seurat_object/integrated_seu.RData"))
} else {
  load(file = paste0(path, "decidua/",type,"_seu.RData"))
}

rowData = as.data.frame(row.names(seu[[assay]]))
colnames(rowData) <- c("feature_symbol")
colData = as.data.frame(colnames(seu[[assay]]))
colnames(colData) = c("Barcode")


if ( origin == "menstrual" ){
  # add seurat cluster number to SingleCellExperiment
  cluster_info = as.data.frame(Idents(seu))
  cluster_info$Barcode = row.names(cluster_info)
  colData = left_join(colData, cluster_info, by = "Barcode")
  colnames(colData)[2] = "cluster"
} else {
  # add celltype to SingleCellExperiment
  colData = read.table("/icgc/dkfzlsdf/analysis/B210/Evelin/E-MTAB-6701_arrayexpress_10x_meta.txt", header = TRUE)
  colData = colData[which(colData$Cell %in% colnames(seu[[assay]])), ]
}


seu_matrix = as.matrix(seu[[assay]]@data)

se = SummarizedExperiment(assays = list(logcounts = matrix(as.numeric(unlist(seu_matrix)),nrow=nrow(seu_matrix))), rowData = rowData, colData = colData)
rm(seu)

if ( origin == "menstrual" ){
  menstrual = as(se, "SingleCellExperiment")
  save(menstrual, file = paste0(path, "sce_RData/",type,"_",origin,".RData"))
} else {
  decidua = as(se, "SingleCellExperiment")
  save(decidua, file = paste0(path, "sce_RData/",type,"_",origin,".RData"))
}