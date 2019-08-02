# First trimester decidua dataset
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6701/files/processed/
decidua_raw = read.table("/icgc/dkfzlsdf/analysis/B210/Evelin/raw_data_10x.txt", header = TRUE, row.names = 1)

# filter out cells that express < 500 genes (UMI?)
# this should not work, but it does
counts = decidua_raw[which(colSums(decidua_raw)>=500)]
rm(decidua_raw)

# alternative for filtering low expression cells that works just as well 
# filter out cells that express < 500 genes
# if value != 0 then value = 1 -> find colsums
# if colSum<500 filter out column
counts = decidua_raw
decidua_raw[decidua_raw!=0]=1
counts = counts[which(colSums(decidua_raw)>=500)]
rm(decidua_raw)

# filter out genes that are expressed in < 3 cells
counts = counts[which(rowSums(counts)>=3)]
write.table(counts, file="/icgc/dkfzlsdf/analysis/B210/Evelin/decidua.txt", col.names = T,quote = F, row.names = T)

# only keep metadata for cells that appear in reduced expression matrix
decidua_meta = read.table("/icgc/dkfzlsdf/analysis/B210/Evelin/E-MTAB-6701_arrayexpress_10x_meta.txt", header = TRUE)
decidua_counts = read.table("/icgc/dkfzlsdf/analysis/B210/Evelin/decidua.txt", header = T, row.names = 1)
decidua_meta = decidua_meta[which(decidua_meta$Cell %in% colnames(decidua_counts)), ]

# choose 5000 of most variable expressed genes
expression_var = apply(decidua_counts, 1, var)
variably_expressed = names(tail(sort(expression_var),5000))
decidua_counts = decidua_counts[variably_expressed, ]

save(decidua_counts, file="/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/decidua_5000_variably_expressed_counts.RData")


