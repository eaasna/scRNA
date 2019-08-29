# samples combined using anchorgenes + SCTransformed count matrix
library(Seurat)
dir_path = "/icgc/dkfzlsdf/analysis/B210/data/mf/"
devtools::install_github(repo = "satijalab/seurat", ref = "develop")
library(future)
options(future.globals.maxSize=891289600)
source("/icgc/dkfzlsdf/analysis/B210/Evelin/git-repo/seurat_utils.R")


norm = "SCT"
load( file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/menstrual_RData/seu_list_",norm,".RData" ))

seu_features <- SelectIntegrationFeatures(object.list = seu_list, nfeatures = 3000)
seu_list <- PrepSCTIntegration(object.list = seu_list, anchor.features = seu_features, verbose = FALSE)

# considering 50 nearest neighbors when filtering anchors <- close to upper limit for smallest sample
anchors <- FindIntegrationAnchors(object.list = seu_list, normalization.method = "SCT", 
                                  anchor.features = seu_features, verbose = FALSE, k.filter = 50)
seu <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE, preserve.order = TRUE )


DefaultAssay(seu) <- "RNA"
nPCA = 20
seu <- NormalizeData(seu, normalization.method = "LogNormalize", assay = "RNA")
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, assay = "RNA")
seu <- ScaleData(seu, features = rownames(seu), vars.to.regress = "percent.mt", assay = "RNA")
seu <- RunPCA(seu, features = VariableFeatures(seu, assay = "RNA"), assay = "RNA", reduction.name = "log.pca")
seu <- RunUMAP(seu, reduction = "log.pca", dims = 1:nPCA, assay = "RNA", reduction.name = "log.umap")

seu <- FindNeighbors(seu, dims = 1:nPCA, reduction = "log.pca", assay = "RNA")
seu <- FindClusters(seu, verbose = FALSE)
seu$log.clusters = Idents(seu)


DefaultAssay(seu) <- "SCT"
seu <- SCTransform(seu, verbose = FALSE, conserve.memory = FALSE, vars.to.regress = "percent.mt", assay = "RNA")
seu <- RunPCA(seu, features = VariableFeatures(seu, assay = "SCT"), assay = "SCT", reduction.name = "pca")
seu <- RunUMAP(seu, reduction = "SCT.pca", dims = 1:nPCA, assay = "SCT", reduction.name = "SCT.umap")

seu <- FindNeighbors(seu, dims = 1:nPCA, reduction = "SCT.pca", assay = "SCT")
seu <- FindClusters(seu, verbose = FALSE)
seu$SCT.clusters = Idents(seu)


save(seu, file = "/icgc/dkfzlsdf/analysis/B210/Evelin/menstrual_RData/integrated_seu.RData")
