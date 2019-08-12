# naively combined samples + log normalized count matrix
library(Seurat)
library(ggplot2)
library(devtools)

dir_path = "/icgc/dkfzlsdf/analysis/B210/data/mf/"
runs = c("cellranger201_count_23156_6_GRCh38", 
         "cellranger201_count_23156_7_GRCh38", 
         "cellranger201_count_23156_8_GRCh38", 
         "cellranger201_count_24192-25608_4839STDY7131581_GRCh38", 
         "cellranger201_count_24192-25608_4839STDY7131582_GRCh38", 
         "cellranger201_count_24192-25608_4839STDY7131583_GRCh38", 
         "cellranger201_count_24192-25608_4839STDY7131584_GRCh38")

seu <- Read10X(paste0(dir_path, runs, "/outs/filtered_gene_bc_matrices/GRCh38/"))
seu <- CreateSeuratObject( seu, min.cells = 3 )

# Storing percentage of mitochondrial genes in object meta data
seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")

#VlnPlot(seu, features = c("percent.mt"))
#ggsave("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/mito_violin.pdf", width = 8, height = 8) + NoLegend()


seu <- NormalizeData(seu, normalization.method = "LogNormalize")
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seu), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seu)
pdf("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/most_variable_genes.pdf", height = 5, width = 5) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + labs(title = "Highly variable genes") +
  theme(plot.title = element_text(size = 18)) + NoLegend()
plot2
dev.off()

VlnPlot(seu, features = c("nFeature_RNA"))
ggsave("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/feature_violin.pdf", width = 8, height = 8) + NoLegend()

#Determining cut off point for low quality (empty) cells
#y <- number of cells (barcode)
#x <- number of expressed genes
n_detected_genes_sorted <- sort(seu$nFeature_RNA, decreasing = TRUE)

#pdf("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/joined_kneeplot_SCT.pdf") 
#plot(seq.int(length(n_detected_genes_sorted)), n_detected_genes_sorted, log = "xy", type = "s", xlab = "barcodes", ylab = "detected genes")
#dev.off()

# n is based on kneeplot
n = 500

# Exclude cells that expressed < n genes
low_exp_cells = names(seu$nFeature_RNA[seu$nFeature_RNA < n])
seu = subset(seu, cells = low_exp_cells, invert = TRUE)

all_genes = rownames(seu)
seu <- ScaleData(seu, features = all_genes, vars.to.regress = "percent.mt")
seu <- RunPCA(seu, features = VariableFeatures(seu))

#finding optimal number of principal components
#ElbowPlot(seu)
#ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/",type,"joined_elbowplot.pdf"), width = 7, height = 5)

# based on elbowplot
nPCA = 30
seu <- RunUMAP(seu, dims = 1:nPCA)
#seu <- RunTSNE(seu, dims = 1:nPCA)

seu <- FindNeighbors(seu, dims = 1:nPCA)
seu <- FindClusters(seu, verbose = FALSE)


save(seu, file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/log_seu.RData"))