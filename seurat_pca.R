library(Seurat)
library(ggplot2)

load(file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/joined.RData"))

# n is based on kneeplot
n = 500

# Exclude cells that expressed < n genes
low_exp_cells = names(seu$nFeature_RNA[seu$nFeature_RNA < n])
seu = subset(seu, cells = low_exp_cells, invert = TRUE)

all_genes = rownames(seu)
seu <- ScaleData(seu, features = all_genes, vars.to.regress = "percent.mt")

seu <- RunPCA(seu, features = VariableFeatures(seu))

#finding optimal number of principal components
ElbowPlot(seu)
ggsave("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/joined_elbowplot.pdf")

save(seu, file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/joined.RData"))
