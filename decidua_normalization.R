library(Seurat)
library(devtools)

path = "/icgc/dkfzlsdf/analysis/B210/Evelin/decidua/"
seu <- Read10X(path)
seu <- CreateSeuratObject( seu, min.features = 500, min.cells = 60 )

seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")
high_mt_cells = names(seu$nFeature_RNA[seu$percent.mt > 20])
seu = subset(seu, cells = high_mt_cells, invert = TRUE)

# seurat slot for cell type annotation
metadata = read.table("/icgc/dkfzlsdf/analysis/B210/Evelin/E-MTAB-6701_arrayexpress_10x_meta.txt", header = TRUE)
metadata = metadata[which(metadata$Cell %in% colnames(seu[["RNA"]])), ]
seu$annotation = metadata$annotation


if (FALSE){
  tmp = seu
  seu <- NormalizeData(seu, normalization.method = "LogNormalize" )
  high_var_genes = VariableFeatures(seu)[1:5000]
  seu = subset(seu, features = high_var_genes )
  seu <- ScaleData(seu, features = rownames(seu) )
  save( seu, file = "/icgc/dkfzlsdf/analysis/B210/Evelin/decidua/log_seu_5000.RData" )
  
  nPCA = 20
  seu <- RunPCA(seu, features = VariableFeatures(seu) )
  #seu <- RunUMAP(seu, dims = 1:nPCA)
  seu <- RunTSNE(seu, dims = 1:nPCA )
  seu <- FindNeighbors(seu, dims = 1:nPCA)
  seu <- FindClusters(seu, verbose = FALSE)
  
  seu = tmp
  rm(tmp)
  
}
# size limit: 2048 * 1024^2 
options(future.globals.maxSize=891289600)
seu <- SCTransform(seu, verbose = FALSE, conserve.memory = TRUE )
high_var_genes = VariableFeatures(seu)[1:5000]
seu = subset(seu, features = high_var_genes )

nPCA = 30
seu <- RunPCA(seu, features = VariableFeatures(seu) )
#seu <- RunUMAP(seu, dims = 1:nPCA)
seu <- RunTSNE(seu, dims = 1:nPCA )

seu <- FindNeighbors(seu, dims = 1:nPCA)
seu <- FindClusters(seu, verbose = FALSE)
save( seu, file = "/icgc/dkfzlsdf/analysis/B210/Evelin/decidua/SCT_seu_5000.RData" )