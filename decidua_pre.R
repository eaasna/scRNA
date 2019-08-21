# First trimester decidua dataset
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6701/files/processed/
decidua.raw = read.table("/icgc/dkfzlsdf/analysis/B210/Evelin/raw_data_10x.txt", header = TRUE, row.names = 1)

dims = dim(decidua.raw)

tense.counts = decidua.raw
decidua.raw[decidua.raw!=0]=1
tense.counts = tense.counts[which(rowSums(decidua.raw)>=60), which(colSums(decidua.raw)>=500)]
rm(decidua.raw)

subset_of_features <- function(nr){
  expression_var = apply(tense.counts, 1, var)
  variably_expressed = names(tail(sort(expression_var), nr))
  tense.counts = tense.counts[variably_expressed, ]
}

sample = sample(64735, 10000)
sample = sort(sample)
write.table(sample, file = "/icgc/dkfzlsdf/analysis/B210/Evelin/decidua/sample", quote = FALSE, row.names = FALSE, col.names = FALSE)
tense.counts = tense.counts[,sample]


write.table(colnames(tense.counts), "/icgc/dkfzlsdf/analysis/B210/Evelin/decidua/barcodes.tsv", sep="\n", row.names = FALSE, quote = FALSE, col.names = FALSE)

# row and column indices of non-zero values
non.zero = as.data.frame(which(tense.counts != 0, arr.ind=TRUE))
colnames(non.zero) = c("row", "col")
library(dplyr)
non.zero = arrange(non.zero, non.zero$row, non.zero$col)

vector.counts = as.vector(t(tense.counts))
vector.counts = vector.counts[which(vector.counts!=0)]

# sparse matrix for pedestrian analysis
#library(Matrix)
#sparse.counts <- sparseMatrix(i = non.zero$row, j = non.zero$col, x = vector.counts, dims = dim(tense.counts), dimnames = dimnames(tense.counts))

# Seurat matrix
# cols: gene i.e row nr, barcode i.e col nr, count
# first non commented row agrees with number of rows in barcode.tsv and genes.tsv
# i.e first non commented row shows total number of genes, barcodes, UMIs 
mtx = data.frame(genes = c(dims[1], non.zero$row), barcodes = c(dims[2], non.zero$col), counts = c(sum(vector.counts), vector.counts))
rm(non.zero, tense.counts)
write.table(mtx, "/icgc/dkfzlsdf/analysis/B210/Evelin/decidua/matrix.mtx", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)