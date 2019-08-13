devtools::install_github("hms-dbmi/conos")
library(conos)
library(Seurat)

norm = "log"
load( file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/seu_list_",norm,".RData" ))

# how many cores on HPC?
con <- Conos$new(seu_list, n.cores=4)
