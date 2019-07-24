df = as.data.frame(seu$orig.ident)
df = as.data.frame(table(df))

titlesize = 20
slices <- df$Freq 
lbls <- df$df
pct <- round(slices/sum(slices)*100)
lbls <- paste0(pct, "%") # add percentage to labels 
#lbls <- paste(lbls, "\n", df$df, sep = "")

pdf("/icgc/dkfzlsdf/analysis/B210/Evelin/plots/pie_afterQC.pdf")
pie(slices,labels = lbls, col=rainbow(length(lbls)), main="Cell origin by individual (after QC)")
dev.off()


df_before_qc = data.frame(df = c(1, 2, 3, 4), Freq = c(737280, 737280, 1474560, 2211840))
df_after_qc_min2000 = data.frame(df = c(1, 2, 3, 4), Freq = c(375, 40, 37, 1029))
df_after_qc_min1000 = data.frame(df=c(1, 2, 3, 4), Freq = c(552, 69, 63, 2265))
df_after_qc_min500 = data.frame(df=c(1, 2, 3, 4), Freq = c(678, 87, 255, 3514))


