library(pROC)
load(file = paste0("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/joined.RData"))


load("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/selPoints.RData")

# lasso slot - boolean variable marks selected points
seu$lasso = 0
seu$lasso[selPoints]=1

# barcode of each cell in selection
selCells = colnames(seu[["RNA"]]@counts)[selPoints]


DimPlot(seu, reduction = "umap", group.by="lasso") + 
  labs(title = paste("Lasso selection")) +
  theme(plot.title = element_text(size = 20))


# row variance for column-sparse matrix
# method from: https://www.zmbh.uni-heidelberg.de/anders/div/sc_pedestrian.html
colVars_spm <- function( spm ) {
  stopifnot( is( spm, "dgCMatrix" ) )
  ans <- sapply( seq.int(spm@Dim[2]), function(j) {
    mean <- sum( spm@x[ (spm@p[j]+1):spm@p[j+1] ] ) / spm@Dim[1]
    sum( ( spm@x[ (spm@p[j]+1):spm@p[j+1] ] - mean )^2 ) +
      mean^2 * ( spm@Dim[1] - ( spm@p[j+1] - spm@p[j] ) ) } ) / ( spm@Dim[1] - 1 )
  names(ans) <- spm@Dimnames[[2]]
  ans
}

rowVars_spm <- function( spm ) {
  colVars_spm( t(spm) )
}

gene_means <- rowMeans( seu[["RNA"]]@counts )
gene_vars <- rowVars_spm( seu[["RNA"]]@counts )

plot( gene_means, gene_vars / gene_means,
      log = "xy", cex = .3, col = adjustcolor("black", alpha=.3), 
      xlab = "mean", ylab = "variance / mean" )

poisson_vmr <- mean( 1 / colSums( seu[["RNA"]]@counts ) )

abline( h = 1:3 * poisson_vmr, col="lightblue" )

informative_genes = names( which( gene_vars / gene_means  >  3 * poisson_vmr ) )
length(informative_genes)