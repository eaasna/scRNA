# naively combined samples + SCTransormed count matrix
library(Seurat)
dir_path = "/icgc/dkfzlsdf/analysis/B210/data/mf/"
library(future)
options(future.globals.maxSize=891289600)
source("/icgc/dkfzlsdf/analysis/B210/Evelin/git-repo/seurat_utils.R")

dir_path = "/icgc/dkfzlsdf/analysis/B210/data/mf/"
runs = c("cellranger201_count_23156_6_GRCh38", 
         "cellranger201_count_23156_7_GRCh38", 
         "cellranger201_count_23156_8_GRCh38", 
         "cellranger201_count_24192-25608_4839STDY7131581_GRCh38", 
         "cellranger201_count_24192-25608_4839STDY7131582_GRCh38", 
         "cellranger201_count_24192-25608_4839STDY7131583_GRCh38", 
         "cellranger201_count_24192-25608_4839STDY7131584_GRCh38")

seu <- Read10X(paste0(dir_path, runs, "/outs/filtered_gene_bc_matrices/GRCh38/"))
seu <- CreateSeuratObject( seu, min.features = 500, min.cells = 3 )

seu <- exclude_high_mt(seu, 20)

seu <- SCTransform(seu, verbose = FALSE, conserve.memory = FALSE, vars.to.regress = "percent.mt")

seu <- RunPCA(seu, features = VariableFeatures(seu, assay = "SCT"), assay = "SCT")
seu <- after_dim_reduc(seu, 30)

save(seu, file = "/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/SCT_seu.RData")