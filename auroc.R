load(file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/joined.RData"))
load("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/goi.RData")

load("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/selPoints.RData")

# lasso slot - boolean variable marks selected points
seu$lasso = 0
seu$lasso[selPoints]=1

# barcode of each cell in selection
selCells = colnames(seu[["RNA"]]@counts)[selPoints]


DimPlot(seu, reduction = "umap", group.by="lasso") + 
  labs(title = paste("Lasso selection")) +
  theme(plot.title = element_text(size = 20))


library(pROC)

informative_genes = VariableFeatures(seu)
selPoints_bool <-seq.int(ncol(seu[["RNA"]]@counts)) %in% selPoints
# from https://www.zmbh.uni-heidelberg.de/anders/div/sc_pedestrian.html
aucs <- sapply( informative_genes, function(g)
  auc( roc( selPoints_bool, seu[["RNA"]]@counts[g,] ) ) )


marker_genes = names( head( sort( aucs, decreasing = TRUE ), 20 ) )

# print name of celltype and which marker genes overlap (if there is overlap)

find_overlap<- function(known_markers){
  if (length(which(marker_genes %in% known_markers)) > 0){
    print(names(known_markers))
    print(marker_genes[which(marker_genes %in% known_markers)])
  }
}

sapply(goi, function(known_markers) find_overlap)
