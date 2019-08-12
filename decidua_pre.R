# First trimester decidua dataset
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6701/files/processed/
decidua.raw = read.table("/icgc/dkfzlsdf/analysis/B210/Evelin/raw_data_10x.txt", header = TRUE, row.names = 1)

# filter out cells that express < 500 genes
# if value != 0 then value = 1 -> find colsums, rowsums
# if colSum<500 filter out column -> all cells express >500 genes
# filter out genes that are expressed in < 3 cells
tense.counts = decidua.raw
decidua.raw[decidua.raw!=0]=1
tense.counts = tense.counts[which(rowSums(decidua.raw)>=3), which(colSums(decidua.raw)>=500)]
rm(decidua.raw)

# creating sparse matrix to save RAM
non.zero = as.data.frame(which(tense.counts != 0, arr.ind=TRUE))
colnames(non.zero) = c("row", "col")
library(dplyr)
non.zero = arrange(non.zero, non.zero$row, non.zero$col)

vector.counts = as.vector(t(tense.counts))
vector.counts = vector.counts[which(vector.counts!=0)]

library(Matrix)
sparse.counts <- sparseMatrix(i = non.zero$row, j = non.zero$col, x = vector.counts, dims = dim(tense.counts), dimnames = dimnames(tense.counts))
rm(non.zero, tense.counts)


save(sparse.counts, "/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/decidua_sparse_matrix.RData")


if (FALSE){
  # only keep metadata for cells that appear in reduced expression matrix
  decidua_meta = read.table("/icgc/dkfzlsdf/analysis/B210/Evelin/E-MTAB-6701_arrayexpress_10x_meta.txt", header = TRUE)
  decidua_meta = decidua_meta[which(decidua_meta$Cell %in% colnames(counts)), ]
  
  # choose nr of most variable expressed genes
  expression_var = apply(counts, 1, var)
  #nr = as.numeric(args[1])
  nr = 5000
  variably_expressed = names(tail(sort(expression_var), nr))
  decidua_counts = counts[variably_expressed, ]
  rm(counts)
  save(decidua_counts, file="/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/decidua_5000_variably_expressed_counts.RData")
}


