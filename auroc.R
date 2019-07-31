library(Seurat)
library(pROC)

load(file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/joined.RData"))
load("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/goi.RData")

load("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/selPoints.RData")

# lasso slot - boolean variable marks selected points
seu$lasso = 0
seu$lasso[selPoints]=1

# barcode of each cell in selection
selCells = colnames(seu[["RNA"]]@counts)[selPoints]


informative_genes = VariableFeatures(seu)
selPoints_bool <-seq.int(ncol(seu[["RNA"]]@counts)) %in% selPoints
# from https://www.zmbh.uni-heidelberg.de/anders/div/sc_pedestrian.html
aucs <- sapply( informative_genes, function(g)
  auc( roc( selPoints_bool, seu[["RNA"]]@counts[g,] ) ) )

marker_genes = names( head( sort( aucs, decreasing = TRUE ), 40 ) )

cell_annotation = ""
gene_annotation = ""
# print name of celltype and which marker genes overlap (if there is overlap)
for (i in 1:length(goi)){
  known_markers = goi[i]
  if (length(which(marker_genes %in% known_markers)) > 0){
    cell_annotation = names(known_markers)
    gene_annotation = paste(marker_genes[which(marker_genes %in% known_markers)])
    print(cell_annotation)
    print(gene_annotation)
  }
}

DimPlot(seu, reduction = "umap", group.by="lasso") + 
  labs(title = paste0(cell_annotation, ": ", gene_annotation)) +
  theme(plot.title = element_text(size = 10))
ggsave("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/lasso_4.pdf")