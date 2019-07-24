library( sleepwalk )
load(file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/joined.RData"))
sleepwalk( seu@reductions$umap@cell.embeddings, seu@reductions$pca@cell.embeddings )
