library(Seurat)
library(ggplot2)
library(plyr)

type = "cca_"
load(file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/",type,"joined.RData"))

# based on elbowplot
nPCA = 10
seu <- FindNeighbors(seu, dims = 1:nPCA)


if (FALSE){
  metadata = read.table("/icgc/dkfzlsdf/analysis/B210/data/mf/metadata.tsv", header = TRUE)
  
  # adding sample identity
  seu$orig.sample = as.factor(substr(colnames(seu), 1, 1))
  seu$orig.sample = revalue(seu$orig.sample, c("A"="1", "C"="1", "T"="1", "G"="1"))
  
  # adding individual identity
  seu$orig.ident = revalue(seu$orig.sample, c("1"="07", "2"="10", "3"="18", "4"="20", "5"="07", "6"="07", "7"="20"))
  
  # adding type identity
  seu$orig.type = revalue(seu$orig.sample, c("1"="f", "2"="f", "3"="f", "4"="f", "5"="f", "6"="c", "7"="c"))
  
  # renaming sample identity based on metadata
  seu$orig.sample = revalue(seu$orig.sample, c("1"="07efM", "2"="10dfM", "3"="18bfM", "4"="20afM", "5"="07dfM", "6"="07dcM", "7"="20acM"))
}


if (TRUE) {
  seu$orig.ident = seu@meta.data[, "sample"]
}


seu <- FindClusters(seu, resolution = 1)

#Look at these methods later:
#BoldTitle()
#FontSize()
titlesize = 20
seu <- RunUMAP(seu, dims = 1:nPCA)
DimPlot( seu, reduction = "umap" ) + 
  labs(title = "UMAP grouped by clustering results") +
  theme(plot.title = element_text(size = titlesize))
ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/",type,"umap.pdf"))

seu <- RunTSNE(seu, dims = 1:nPCA)
DimPlot(seu, reduction = "tsne") + 
  labs(title = "tSNE grouped by clustering results") +
  theme(plot.title = element_text(size = titlesize))
ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/",type,"tsne.pdf"))


reds = c("tsne", "umap")
labels = c("tSNE", "UMAP")
for (i in range(1,2)){
#if (FALSE){
  #DimPlot(seu, reduction = reds[i], group.by="orig.sample") + 
  #  labs(title = paste(labels[i],"grouped by batch")) +
  #  theme(plot.title = element_text(size = titlesize))
  #ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/", type, reds[i], "_by_sample.pdf"))
  
  
  DimPlot(seu, reduction = reds[i], group.by="orig.ident") + 
    labs(title = paste(labels[i],"grouped by individual")) +
    theme(plot.title = element_text(size = titlesize))
  ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/", type, reds[i], "_by_individual.pdf"))
  
  
  #DimPlot(seu, reduction = reds[i], group.by="orig.type") + 
  #  labs(title = paste(labels[i],"grouped by type")) +
  #  theme(plot.title = element_text(size = titlesize))
  #ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/", type, reds[i],"_by_type.pdf"))
}

save(seu, file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/",type,"joined.RData"))
