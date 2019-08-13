# samples combined using anchorgenes + SCTransformed count matrix
library(Seurat)
dir_path = "/icgc/dkfzlsdf/analysis/B210/data/mf/"
devtools::install_github(repo = "satijalab/seurat", ref = "develop")
library(future)
options(future.globals.maxSize=891289600)

load( file = "/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/seu_list_sct.RData" )

seu_features <- SelectIntegrationFeatures(object.list = seu_list, nfeatures = 5000)
seu_list <- PrepSCTIntegration(object.list = seu_list, anchor.features = seu_features, verbose = FALSE)

# considering 50 nearest neighbors when filtering anchors <- close to upper limit for smallest sample
anchors <- FindIntegrationAnchors(object.list = seu_list, normalization.method = "SCT", 
                                  anchor.features = seu_features, verbose = FALSE, k.filter = 50)
seu <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
seu <- RunPCA(seu, features = VariableFeatures(seu, assay = "integrated"), assay = "integrated")

nPCA = 30
seu <- RunUMAP(seu, dims = 1:nPCA)
#seu <- RunTSNE(seu, dims = 1:nPCA)

seu <- FindNeighbors(seu, dims = 1:nPCA)
seu <- FindClusters(seu, verbose = FALSE)

save(seu, file = "/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/integrated_seu.RData")