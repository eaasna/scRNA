library(Seurat)
dir_path = "/icgc/dkfzlsdf/analysis/B210/data/mf/"

# manually combine samples for each donor
seven = Read10X( paste0( dir_path, c("cellranger201_count_23156_6_GRCh38", "cellranger201_count_24192-25608_4839STDY7131582_GRCh38", "cellranger201_count_24192-25608_4839STDY7131583_GRCh38"), "/outs/filtered_gene_bc_matrices/GRCh38/"))
ten = Read10X( paste0( dir_path, c("cellranger201_count_23156_7_GRCh38"), "/outs/filtered_gene_bc_matrices/GRCh38/"))
eighteen = Read10X( paste0( dir_path, c("cellranger201_count_23156_8_GRCh38"), "/outs/filtered_gene_bc_matrices/GRCh38/"))
twenty = Read10X( paste0( dir_path, c("cellranger201_count_24192-25608_4839STDY7131581_GRCh38", "cellranger201_count_24192-25608_4839STDY7131584_GRCh38"), "/outs/filtered_gene_bc_matrices/GRCh38/"))

seu_list = c(seven, ten, eighteen, twenty)
seu_after_processing = c()
for(seu in seu_list){
  seu <- CreateSeuratObject( seu, min.cells = 3, min.features = 500)
  seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")
  high_mt_cells = names(seu$nFeature_RNA[seu$percent.mt > 20])
  seu = subset(seu, cells = high_mt_cells, invert = TRUE)
  seu <- NormalizeData(seu, normalization.method = "LogNormalize")
  all_genes = rownames(seu)
  seu <- ScaleData(seu, features = all_genes, vars.to.regress = "percent.mt")
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu_after_processing = c(seu_after_processing, seu)
}

rm(seu_list)

seven = seu_after_processing[[1]]
ten = seu_after_processing[[2]]
eighteen = seu_after_processing[[3]]
twenty = seu_after_processing[[4]]

seven@meta.data[,"sample"] <- "07"
ten@meta.data[,"sample"] <- "10"
eighteen@meta.data[,"sample"] <- "18"
twenty@meta.data[,"sample"] <- "20"

rm(seu_after_processing)

# use CCA to combine all cells between donors
seuA = RunCCA(object1 = seven, object2 = ten)
seuB = RunCCA(object1 = eighteen, object2 = twenty)
seuA = FindVariableFeatures(seuA, selection.method = "vst", nfeatures = 2000)
seuB = FindVariableFeatures(seuB, selection.method = "vst", nfeatures = 2000)
seu = RunCCA(object1 = seuA, object2 = seuB)

# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = seu, reduction.use = "cca", group.by = "sample", pt.size = 0.5, 
              do.return = TRUE)
p2 <- VlnPlot(object = seu, group.by = "sample", do.return = TRUE, features = "nFeature_RNA")
library(cowplot)
plot_grid(p1, p2)


anchors <- FindIntegrationAnchors(object.list = list(seven, ten, eighteen, twenty), dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)

