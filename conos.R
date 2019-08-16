devtools::install_github("hms-dbmi/conos")
library(conos)
library(Seurat)

norm = "SCT"
load( file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/seu_list_",norm,".RData" ))
load( file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/decidua/",norm,"_seu_5000.RData"))
seu_list[[length(seu_list)+1]] = seu
rm(seu)


names(seu_list) = c("23156_6", "24192-25608_4839STDY7131582", "24192-25608_4839STDY7131583", "23156_7", 
                    "23156_8", "24192-25608_4839STDY7131581", "24192-25608_4839STDY7131584", "decidua")


seu_list[[6]] <- NULL
seu_list[[5]] <- NULL
# saving seu_list that has 2 samples removed
# TSNE error: perplexity too large for the number of samples
save(seu_list, file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/seu_list_",norm,"_test.RData"))

# how many cores on HPC?
con <- Conos$new(seu_list, n.cores=4)

#clusters based on each sample separately
pdf(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/con_clusters.pdf"))
con$plotPanel(clustering="multilevel", use.local.clusters=T, title.size=6)
dev.off()

con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30, n.odgenes=2000, 
               matching.method='mNN', metric='angular', 
               score.component.variance=TRUE, verbose=TRUE)

con$findCommunities(method=leiden.community, resolution=1)

#corresponding clusters between samples
pdf(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/con_united_clusters.pdf"))
con$plotPanel(font.size=4)
dev.off()

cellannot <- seu$annotation



# propagating labels from reference
new.label.probabilities <- con$propagateLabels(labels = cellannot, verbose=T)

pdf(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/uncertainty.pdf"))
con$plotGraph(colors=(1 - apply(new.label.probabilities, 1, max)), show.legend=T, legend.title="Uncertainty", legend.pos=c(1, 0))
dev.off()

new.annot <- setNames(colnames(new.label.probabilities)[apply(new.label.probabilities,1,which.max)], rownames(new.label.probabilities))

#clusters labelled based on reference
pdf(paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/con_labelled_clusters.pdf"))
con$plotPanel(groups = new.annot)
dev.off()


write.table(as.data.frame(new.annot), file = "/icgc/dkfzlsdf/analysis/B210/Evelin/new_annot_sct.txt", col.names = F, quote = F, row.names = T, append = F)

