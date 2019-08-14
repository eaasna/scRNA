library(Seurat)

# return Seurat object without cells expressing more than percentage mitochondrial genes
exclude_high_mt <- function(seu, percentage){
  seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")
  high_mt_cells = names(seu$nFeature_RNA[seu$percent.mt > percentage])
  seu = subset(seu, cells = high_mt_cells, invert = TRUE)
  return(seu)
}

# return Seurat object with UMAP embedding, neighbors, clusters
after_dim_reduc <- function(seu, nPCA){
  seu <- RunUMAP(seu, dims = 1:nPCA)
  #seu <- RunTSNE(seu, dims = 1:nPCA)
  
  seu <- FindNeighbors(seu, dims = 1:nPCA)
  seu <- FindClusters(seu, verbose = FALSE)
}


qc_plots <- function(seu, suffix){
  VlnPlot(seu, features = c("percent.mt"))
  ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/",suffix,"mito_violin.pdf"), width = 8, height = 8) + NoLegend()
  
  VlnPlot(seu, features = c("nFeature_RNA"))
  ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/",suffix,"feature_violin.pdf"), width = 8, height = 8) + NoLegend()
  
  #Determining cut off point for low quality (empty) cells
  #y <- number of cells (barcode)
  #x <- number of expressed genes
  n_detected_genes_sorted <- sort(seu$nFeature_RNA, decreasing = TRUE)
  pdf(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/",suffix,"_kneeplot_SCT.pdf")) 
  plot(seq.int(length(n_detected_genes_sorted)), n_detected_genes_sorted, log = "xy", type = "s", xlab = "barcodes", ylab = "detected genes")
  dev.off()
  
  #finding optimal number of principal components
  ElbowPlot(seu)
  ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/",suffix,"joined_elbowplot.pdf"), width = 7, height = 5)
  
}