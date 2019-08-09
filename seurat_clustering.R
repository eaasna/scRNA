library(Seurat)
library(ggplot2)
library(plyr)

type = "sct_cca_"
load(file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/",type,"joined.RData"))



# based on elbowplot
nPCA = 30
seu <- RunUMAP(seu, dims = 1:nPCA)
seu <- RunTSNE(seu, dims = 1:nPCA)

seu <- FindNeighbors(seu, dims = 1:nPCA)
seu <- FindClusters(seu, verbose = FALSE)

seu$orig.ident = seu@meta.data[, "sample"]

#Look at these methods later:
#BoldTitle()
#FontSize()
titlesize = 20

DimPlot( seu, reduction = "umap" ) + 
  labs(title = "UMAP grouped by clustering results") +
  theme(plot.title = element_text(size = titlesize))
ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/",type,"umap.pdf"), height = 5, width = 5)


DimPlot(seu, reduction = "tsne" ) + 
  labs(title = "tSNE grouped by clustering results") +
  theme(plot.title = element_text(size = titlesize)) 
ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/",type,"tsne.pdf"), height = 5, width = 5)


reds = c("tsne", "umap")
labels = c("tSNE", "UMAP")
for (i in range(1,2)){
#if (FALSE){
  DimPlot(seu, reduction = reds[i], group.by="orig.ident" ) + 
    labs(title = paste(labels[i],"grouped by individual")) +
    theme(plot.title = element_text(size = titlesize)) 
  ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/", type, reds[i], "_by_individual.pdf"), height = 5, width = 5)
}

save(seu, file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/",type,"joined.RData"))
