#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(scmap)
library(SingleCellExperiment)

assay = args[1]

load(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/sce_RData/",assay,"_decidua.RData"))
load(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/sce_RData/",assay,"_menstrual.RData"))


decidua <- selectFeatures(decidua, suppress_plot = TRUE, 100)
decidua <- indexCluster(decidua, cluster_col = "annotation")

library(pheatmap)
pdf(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/",assay,"50_variable_decidua.pdf"))
pheatmap(metadata(decidua)$scmap_cluster_index, show_rownames = FALSE)
dev.off()


scmapCluster_results <- scmapCluster(
  projection = menstrual, 
  index_list = list(
    decidua = metadata(decidua)$scmap_cluster_index
  )
)

save(scmapCluster_results, paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/",assay,"_cluster_results.RData"))

df = as.data.frame(table(colData(decidua)$annotation))
ggplot(df, aes(x = reorder(Var1, -Freq), y = Freq)) + geom_col() + theme_bw() + theme(axis.text.x = element_text(angle = 45)) + labs(title = "decidua reference celltypes")
ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/",assay,"_reference_celltypes.pdf"))


df = as.data.frame(table(as.data.frame(scmapCluster_results)$decidua))
ggplot(df, aes(x = reorder(Var1, -Freq), y = Freq)) + geom_col() + theme_bw() + theme(axis.text.x = element_text(angle = 45)) + labs(title = "celltypes based on all genes")
ggsave(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/",assay,"_celltypes_all_genes.pdf"))


