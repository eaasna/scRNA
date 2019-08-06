# First trimester decidua dataset
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6701/files/processed/
decidua_raw = read.table("/icgc/dkfzlsdf/analysis/B210/Evelin/raw_data_10x.txt", header = TRUE, row.names = 1)

print(dim(decidua_raw))
# filter out cells that express < 500 genes
# if value != 0 then value = 1 -> find colsums, rowsums
# if colSum<500 filter out column -> all cells express >500 genes
# filter out genes that are expressed in < 3 cells
counts = decidua_raw
decidua_raw[decidua_raw!=0]=1
counts = counts[which(rowSums(decidua_raw)>=3), which(colSums(decidua_raw)>=500)]

print(dim(counts))
rm(decidua_raw)


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


