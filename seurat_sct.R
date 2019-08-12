# naively combined samples + SCTransormed count matrix
library(Seurat)
dir_path = "/icgc/dkfzlsdf/analysis/B210/data/mf/"
library(future)
options(future.globals.maxSize=891289600)

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

seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")
high_mt_cells = names(seu$nFeature_RNA[seu$percent.mt > 20])
seu = subset(seu, cells = high_mt_cells, invert = TRUE)

seu <- SCTransform(seu, verbose = FALSE, conserve.memory = FALSE, vars.to.regress = "percent.mt")

seu <- RunPCA(seu, features = VariableFeatures(seu, assay = "SCT"), assay = "SCT")

nPCA = 30
seu <- RunUMAP(seu, dims = 1:nPCA)
#seu <- RunTSNE(seu, dims = 1:nPCA)

seu <- FindNeighbors(seu, dims = 1:nPCA)
seu <- FindClusters(seu, verbose = FALSE)

save(seu, file = "/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/SCT_seu.RData")