library(Seurat)
library(dplyr)
library(ggplot2)

load(file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/joined.RData"))
load("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/goi.RData")

all_markers = unlist(goi, use.names = F)
variable_genes = VariableFeatures(seu)
variable_markers = all_markers[which(all_markers %in% variable_genes)]
celltypes = names(goi)

plotted_celltypes = c()


celltype = "up_stroma_vs_gland"
load("/icgc/dkfzlsdf/analysis/B210/Evelin/plotted_celltypes")
plotted_celltypes[length(plotted_celltypes)+1] = celltype
celltypes = celltypes[which(!(celltypes %in% plotted_celltypes))]
save(plotted_celltypes, file="/icgc/dkfzlsdf/analysis/B210/Evelin/plotted_celltypes")
markers_to_plot = goi[[celltype]][which(goi[[celltype]] %in% variable_markers)]
#markers_to_plot = goi[[celltype]]
print(length(markers_to_plot))

col = 3
h = 15
w = 15
VlnPlot(seu, features = markers_to_plot, ncol = col, pt.size = 0.1) + 
  labs(title = celltype) +
  theme(plot.title = element_text(size = 15))
ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/vln_",celltype,".pdf"), height = h, width = w)

FeaturePlot(seu, features = markers_to_plot, reduction = "umap", ncol = col, pt.size = 0.2) + 
  labs(title = celltype) +
  theme(plot.title = element_text(size = 15))
ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/umap_",celltype,".pdf"), height = h, width = w)



