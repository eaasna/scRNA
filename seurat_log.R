# naively combined samples + log normalized count matrix
library(Seurat)
library(ggplot2)
library(devtools)
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
seu <- CreateSeuratObject( seu, min.cells = 3, min.features = 500 )

# Storing percentage of mitochondrial genes in object meta data
seu <- exclude_high_mt(seu, 20)

seu <- NormalizeData(seu, normalization.method = "LogNormalize")
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

# Exclude cells that expressed < n genes. Based on kneeplot
#n = 500
#low_exp_cells = names(seu$nFeature_RNA[seu$nFeature_RNA < n])
#seu = subset(seu, cells = low_exp_cells, invert = TRUE)

all_genes = rownames(seu)
seu <- ScaleData(seu, features = rownames(seu), vars.to.regress = "percent.mt")
seu <- RunPCA(seu, features = VariableFeatures(seu))

seu <- after_dim_reduc(seu, 30)

save(seu, file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/log_seu.RData"))