#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(Seurat)
library(future)
source("/icgc/dkfzlsdf/analysis/B210/Evelin/git-repo/seurat_utils.R")

dir_path = "/icgc/dkfzlsdf/analysis/B210/data/mf/"
suffix = "/outs/filtered_gene_bc_matrices/GRCh38/"
seven1 = Read10X( paste0( dir_path, c("cellranger201_count_23156_6_GRCh38"), suffix))
seven2 = Read10X(paste0( dir_path, c("cellranger201_count_24192-25608_4839STDY7131582_GRCh38"), suffix))
seven3 = Read10X(paste0( dir_path, c("cellranger201_count_24192-25608_4839STDY7131583_GRCh38"), suffix))
ten = Read10X( paste0( dir_path, c("cellranger201_count_23156_7_GRCh38"), suffix))
eighteen = Read10X( paste0( dir_path, c("cellranger201_count_23156_8_GRCh38"), suffix))
twenty1 = Read10X( paste0( dir_path, c("cellranger201_count_24192-25608_4839STDY7131581_GRCh38"), suffix))
twenty2 = Read10X( paste0( dir_path, c("cellranger201_count_24192-25608_4839STDY7131584_GRCh38"), suffix))

seu_list = c(seven1, seven2, seven3, ten, eighteen, twenty1, twenty2)
annot_list = c("07", "07", "07", "10", "18", "20", "20")

norm = args[1]

for(i in 1:length(seu_list)){
  
  seu_list[[i]] <- CreateSeuratObject( seu_list[[i]], min.features = 500 )
  seu_list[[i]]@meta.data[,"sample"] <- annot_list[i]
  
  seu_list[[i]] <- exclude_high_mt(seu_list[[i]], 20)
  
  if (norm == "SCT"){
    options(future.globals.maxSize=891289600)
    seu_list[[i]] <- SCTransform(seu_list[[i]], verbose = FALSE, conserve.memory = FALSE, vars.to.regress = "percent.mt", new.assay.name = "SCT.single")

  } else {
    seu_list[[i]] <- NormalizeData(seu_list[[i]], normalization.method = "LogNormalize")
    seu_list[[i]] <- FindVariableFeatures(seu_list[[i]], selection.method = "vst", nfeatures = 2000, assay = "RNA")
    all_genes <- rownames(seu_list[[i]])
    seu_list[[i]] <- ScaleData(seu_list[[i]], features = all_genes, vars.to.regress = "percent.mt")

  }
  
  nPCA = 30
  seu_list[[i]] <- RunPCA(seu_list[[i]], features = VariableFeatures(seu_list[[i]]) )
  #seu_list[[i]] <- RunUMAP(seu_list[[i]], dims = 1:nPCA)
  
  if (!(i %in% c(5, 6))){
    seu_list[[i]] <- RunTSNE(seu_list[[i]])  
  }
  
  seu_list[[i]] <- FindNeighbors(seu_list[[i]])
  seu_list[[i]] <- FindClusters(seu_list[[i]], verbose = FALSE)
  
}

load( file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/decidua/",norm,"_seu.RData" ))
seu_list[[8]] = seu
seu_list[[8]]@meta.data[, "sample"] <- "decidua"
rm(seu)

save( seu_list, file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/seu_list_",norm,".RData" ))