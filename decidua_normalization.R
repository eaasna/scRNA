library(Seurat)
library(devtools)

path = "/icgc/dkfzlsdf/analysis/B210/Evelin/decidua/"
seu <- Read10X(path)
seu <- CreateSeuratObject( seu, min.features = 500, min.cells = 3 )


seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")
high_mt_cells = names(seu$nFeature_RNA[seu$percent.mt > 20])
seu = subset(seu, cells = high_mt_cells, invert = TRUE)

# seurat slot for cell type annotation
metadata = read.table("/icgc/dkfzlsdf/analysis/B210/Evelin/E-MTAB-6701_arrayexpress_10x_meta.txt", header = TRUE)
metadata = metadata[which(metadata$Cell %in% colnames(seu[["RNA"]])), ]
seu$annotation = metadata$annotation


tmp = seu
seu <- NormalizeData(seu, normalization.method = "LogNormalize")
all_genes <- rownames(seu)
seu <- ScaleData(seu, features = all_genes, vars.to.regress = "percent.mt")
save( seu, file = "/icgc/dkfzlsdf/analysis/B210/Evelin/decidua/log_seu.RData" )

seu = tmp
rm(tmp)
# size limit: 2048 * 1024^2 
options(future.globals.maxSize=891289600)
seu <- SCTransform(seu, verbose = FALSE, conserve.memory = TRUE, vars.to.regress = "percent.mt")
save( seu, file = "/icgc/dkfzlsdf/analysis/B210/Evelin/decidua/SCT_seu.RData" )

