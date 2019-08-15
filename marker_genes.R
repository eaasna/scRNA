library(Seurat)
library(dplyr)
library(ggplot2)
source("/icgc/dkfzlsdf/analysis/B210/Evelin/git-repo/marker_library.R")


load(file = "/icgc/dkfzlsdf/analysis/B210/Evelin/integrated_seu.RData")
#load("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/goi.RData")
goi = marker.list


#preparing clustering based on assay
assay = "SCT"
type = "SCT"
nPCA = 30

DefaultAssay(seu) <- assay
seu <- FindNeighbors(seu, dims = 1:nPCA, reduction = paste0(type, ".pca"), assay = assay)
seu <- FindClusters(seu, verbose = FALSE)


titlesize = 15
DimPlot( seu, reduction = "log.umap" ) + 
  labs(title = type) +
  theme(plot.title = element_text(size = titlesize))
ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/",type,"_umap.pdf"), height = 5, width = 5)


df = data.frame(col = c(1, 2, 3, 2, 3, 3, 3), h = c(5, 5, 5, 10, 10, 10, 15), w = c(5, 10, 15, 10, 15, 15, 15))

all_markers = unlist(goi, use.names = F)
variable_genes = VariableFeatures(seu)
variable_markers = all_markers[which(all_markers %in% variable_genes)]
celltypes = names(goi)

for (celltype in names(marker.list)){
  #markers_to_plot = goi[[celltype]][which(marker.list[[celltype]] %in% variable_markers)]
  markers_to_plot = marker.list[[celltype]][which(marker.list[[celltype]] %in% row.names(seu[[assay]]@data))]
  # number of marker genes that have UMI counts in this assay
  len = length(markers_to_plot)
  if (len>0 & len<8){
    col = df[len,"col"]
    h = df[len,"h"]
    w = df[len,"w"]
    
    if (len==1) {title = markers_to_plot[1]} else {title = celltype}
    VlnPlot(seu, features = markers_to_plot, ncol = col, pt.size = 0.1) + 
      labs(title = title) +
      theme(plot.title = element_text(size = 15))
    ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/",type,"/vln_",celltype,"_",type,".pdf"), height = h, width = w)
    
    FeaturePlot(seu, features = markers_to_plot, reduction = "log.umap", ncol = col, pt.size = 0.2) + 
      labs(title = title) +
      theme(plot.title = element_text(size = 15))
    ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/",type,"/umap_",celltype,"_",type,".pdf"), height = h, width = w)
  }
}
