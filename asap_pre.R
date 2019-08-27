library(Seurat)
source("/icgc/dkfzlsdf/analysis/B210/Evelin/git-repo/seurat_utils.R")

norm = "SCT"
load(file = "/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/integrated_seu.RData")

if (norm == "log"){assay = "RNA"}
if (norm == "SCT"){assay = "SCT"}

write.table(as.matrix(seu[[assay]]@counts),
            file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/",norm,"_tense_matrix.txt"), 
            col.names = TRUE, quote = FALSE, row.names = TRUE, sep = "\t")


# a file with two columns 
# first column: list of barcodes
# second column, sample ID of that cell
write.table(data.frame(barcode = colnames(seu[[assay]]), sample = as.vector(seu$sample)),
            file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/",norm,"_batch_file.txt"),
            col.names = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")