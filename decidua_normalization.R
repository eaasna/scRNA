library(Seurat)

path = "/icgc/dkfzlsdf/analysis/B210/Evelin/decidua/"
seu <- Read10X(path)
seu <- CreateSeuratObject( seu )

seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")
high_mt_cells = names(seu$nFeature_RNA[seu$percent.mt > 20])
seu = subset(seu, cells = high_mt_cells, invert = TRUE)


seu_sct <- SCTransform(seu, verbose = FALSE, conserve.memory = FALSE, vars.to.regress = "percent.mt")
save( seu_sct, file = "/icgc/dkfzlsdf/analysis/B210/Evelin/decidua/sct_seu.RData" )


seu_log <- NormalizeData(seu, normalization.method = "LogNormalize")
all_genes <- rownames(seu_log)
seu_log <- ScaleData(seu_log, features = all_genes, vars.to.regress = "percent.mt")
save( seu_log, file = "/icgc/dkfzlsdf/analysis/B210/Evelin/decidua/log_seu.RData" )

