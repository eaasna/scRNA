# samples combined using anchorgenes + SCTransformed count matrix
library(Seurat)
dir_path = "/icgc/dkfzlsdf/analysis/B210/data/mf/"
devtools::install_github(repo = "satijalab/seurat", ref = "develop")
library(future)
options(future.globals.maxSize=891289600)

suffix = "/outs/filtered_gene_bc_matrices/GRCh38/"
seven1 = Read10X( paste0( dir_path, c("cellranger201_count_23156_6_GRCh38"), suffix))
seven2 = Read10X(paste0( dir_path, c("cellranger201_count_24192-25608_4839STDY7131582_GRCh38"), suffix))
seven3 = Read10X(paste0( dir_path, c("cellranger201_count_24192-25608_4839STDY7131583_GRCh38"), suffix))
ten = Read10X( paste0( dir_path, c("cellranger201_count_23156_7_GRCh38"), suffix))
eighteen = Read10X( paste0( dir_path, c("cellranger201_count_23156_8_GRCh38"), suffix))
twenty1 = Read10X( paste0( dir_path, c("cellranger201_count_24192-25608_4839STDY7131581_GRCh38"), suffix))
twenty2 = Read10X( paste0( dir_path, c("cellranger201_count_24192-25608_4839STDY7131584_GRCh38"), suffix))


seu_list = c(seven1, seven2, seven3, ten, eighteen, twenty1, twenty2)
annot_list = c("07", "07", "07", "10", "18", "20", "20")
for(i in 1:length(seu_list)){
  
  seu_list[[i]] <- CreateSeuratObject( seu_list[[i]], min.features = 500 )
  seu_list[[i]]@meta.data[,"sample"] <- annot_list[i]
  
  seu_list[[i]] <- PercentageFeatureSet(seu_list[[i]], pattern = "^MT-", col.name = "percent.mt")
  high_mt_cells = names(seu_list[[i]]$nFeature_RNA[seu_list[[i]]$percent.mt > 20])
  seu_list[[i]] = subset(seu_list[[i]], cells = high_mt_cells, invert = TRUE)
  
  # SCTransform replaces NormalizeData, FindVariableFeatures, ScaleData
  # DO NOT run ScaleData after SCTransform
  seu_list[[i]] <- SCTransform(seu_list[[i]], verbose = FALSE, conserve.memory = FALSE, vars.to.regress = "percent.mt")
}

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