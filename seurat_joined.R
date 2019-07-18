library(Seurat)
library(dplyr)
library(ggplot2)

dir_path = "/icgc/dkfzlsdf/analysis/B210/data/mf/"
runs = c("cellranger201_count_23156_6_GRCh38", "cellranger201_count_23156_7_GRCh38", "cellranger201_count_23156_8_GRCh38", "cellranger201_count_24192-25608_4839STDY7131581_GRCh38", "cellranger201_count_24192-25608_4839STDY7131582_GRCh38", "cellranger201_count_24192-25608_4839STDY7131583_GRCh38", "cellranger201_count_24192-25608_4839STDY7131584_GRCh38")

seu <- Read10X(paste0(dir_path, runs, "/outs/raw_gene_bc_matrices/GRCh38/"))
seu <- CreateSeuratObject( seu, min.cells = 3 )
seu <- NormalizeData(seu)


#Determining cut off point for low quality (empty) cells
#y <- number of cells (barcode)
#x <- number of expressed genes
n_detected_genes_sorted <- sort(seu$nFeature_RNA, decreasing = TRUE)

#pdf("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/joined_kneeplot.pdf") 
#plot(seq.int(length(n_detected_genes_sorted)), n_detected_genes_sorted, log = "xy", type = "s", xlab = "barcodes", ylab = "detected genes")
#dev.off()

n = 2000
low_exp_cells = names(seu$nFeature_RNA[seu$nFeature_RNA < n])
seu = subset(seu, cells = low_exp_cells, invert = TRUE)

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
all_genes = rownames(seu)
seu <- ScaleData(seu, features = all_genes, vars.to.regress = NULL)
seu <- RunPCA(seu, features = VariableFeatures(seu))

#finding optimal number of principal components
#ElbowPlot(seu)
#ggsave("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/joined_elbowplot.pdf")

nPCA = 15
seu <- FindNeighbors(seu, dims = 1:nPCA)
seu <- FindClusters(seu, resolution = 1)

seu <- RunUMAP(seu, dims = 1:nPCA)
DimPlot( seu, reduction = "umap" )
ggsave("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/umap_joined.pdf")

seu <- RunTSNE(seu, dims = 1:7)
DimPlot(seu, reduction = "tsne")
ggsave("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/tsne_joined.pdf")

save(seu, file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/joined.RData"))