library(Seurat)
library(ggplot2)

load(file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/joined.RData"))

# based on elbowplot
nPCA = 15

seu <- FindNeighbors(seu, dims = 1:nPCA)

seu <- RunUMAP(seu, dims = 1:nPCA)
DimPlot( seu, reduction = "umap" )
ggsave("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/umap_joined.pdf")

seu <- RunTSNE(seu, dims = 1:nPCA)
DimPlot(seu, label = TRUE, reduction = "tsne") + NoLegend()
ggsave("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/tsne_joined.pdf")

seu <- FindClusters(seu, resolution = 1)

save(seu, file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/joined.RData"))