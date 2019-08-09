library(Seurat)
dir_path = "/icgc/dkfzlsdf/analysis/B210/data/mf/"
devtools::install_github(repo = "satijalab/seurat", ref = "develop")
library(future)
options(future.globals.maxSize=891289600)

# manually combine samples for each donor
seven1 = Read10X( paste0( dir_path, c("cellranger201_count_23156_6_GRCh38"), "/outs/filtered_gene_bc_matrices/GRCh38/"))
seven2 = Read10X(paste0( dir_path, c("cellranger201_count_24192-25608_4839STDY7131582_GRCh38"), "/outs/filtered_gene_bc_matrices/GRCh38/"))
seven3 = Read10X(paste0( dir_path, c("cellranger201_count_24192-25608_4839STDY7131583_GRCh38"), "/outs/filtered_gene_bc_matrices/GRCh38/"))

ten = Read10X( paste0( dir_path, c("cellranger201_count_23156_7_GRCh38"), "/outs/filtered_gene_bc_matrices/GRCh38/"))

eighteen = Read10X( paste0( dir_path, c("cellranger201_count_23156_8_GRCh38"), "/outs/filtered_gene_bc_matrices/GRCh38/"))

twenty1 = Read10X( paste0( dir_path, c("cellranger201_count_24192-25608_4839STDY7131581_GRCh38"), "/outs/filtered_gene_bc_matrices/GRCh38/"))
twenty2 = Read10X( paste0( dir_path, c("cellranger201_count_24192-25608_4839STDY7131584_GRCh38"), "/outs/filtered_gene_bc_matrices/GRCh38/"))


seu_list = c(seven1, seven2, seven3, ten, eighteen, twenty1, twenty2)
annot_list = c("07", "07", "07", "10", "18", "20", "20")
for(i in 1:length(seu_list)){
  
# should actually not be pre-filtering based on no of min.cells expressing each gene
# min.features = 500, min.cells = NA 
# PrepSCTIntegration ERROR: exceeded maximum allowed size of 500 MiB
# workaround: discarding more data before integration
  seu_list[[i]] <- CreateSeuratObject( seu_list[[i]], min.features = 500 )
  seu_list[[i]]@meta.data[,"sample"] <- annot_list[i]
  
  seu_list[[i]] <- PercentageFeatureSet(seu_list[[i]], pattern = "^MT-", col.name = "percent.mt")
  high_mt_cells = names(seu_list[[i]]$nFeature_RNA[seu_list[[i]]$percent.mt > 20])
  seu_list[[i]] = subset(seu_list[[i]], cells = high_mt_cells, invert = TRUE)
  
  # SCTransform replaces NormalizeData, FindVariableFeatures, ScaleData
  # DO NOT run ScaleData after SCTransform
  seu_list[[i]] <- SCTransform(seu_list[[i]], verbose = FALSE, conserve.memory = FALSE, vars.to.regress = "percent.mt")
}

seu_features <- SelectIntegrationFeatures(object.list = seu_list, nfeatures = 3000)
seu_list <- PrepSCTIntegration(object.list = seu_list, anchor.features = seu_features, verbose = FALSE)


# considering 80 nearest neighbors when filtering anchors <- close to upper limit for smallest sample
anchors <- FindIntegrationAnchors(object.list = seu_list, normalization.method = "SCT", anchor.features = seu_features, verbose = FALSE, k.filter = 80)
seu <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)


seu <- RunPCA(seu, features = VariableFeatures(seu))

save(seu, file = "/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/sct_cca_joined.RData")

# Follow this by seurat_clustering.R