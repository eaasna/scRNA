#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


#norm = args[1]
norm = "SCT"


library(SingleCellExperiment)
library(SummarizedExperiment)
library(Seurat)
library(dplyr)
source("/icgc/dkfzlsdf/analysis/B210/Evelin/git-repo/seurat_utils.R")

path = "/icgc/dkfzlsdf/analysis/B210/Evelin/"


if (norm == "log") {assay = "RNA"}
if (norm == "SCT") {assay = "SCT"}


load(file = paste0(path, "decidua/",norm,"_seu.RData"))

rowData = data.frame(feature_symbol = row.names(seu[[assay]]))
colData = data.frame(Barcode = colnames(seu[[assay]]))

# add celltype to SingleCellExperiment
colData = read.table("/icgc/dkfzlsdf/analysis/B210/Evelin/E-MTAB-6701_arrayexpress_10x_meta.txt", header = TRUE)
colData = colData[which(colData$Cell %in% colnames(seu[[assay]])), ]


assays = return_assay(seu, assay)

se = SummarizedExperiment(assays = list(counts = assays[[1]], logcounts = assays[[2]]), rowData = rowData, colData = colData)
rm(seu)

decidua = as(se, "SingleCellExperiment")
save(decidua, file = paste0(path, "sce_RData/", norm,"_decidua.RData"))
