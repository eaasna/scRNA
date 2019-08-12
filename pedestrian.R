library( pROC )
load("/icgc/dkfzlsdf/analysis/B210/Evelin/seurat_object/decidua_sparse_matrix.RData")

# normalizing for sequnecing depth
nrm_counts <- t( t(counts) / colSums(counts) )

# finding highly variable genes
# https://www.zmbh.uni-heidelberg.de/anders/div/sc_pedestrian.html
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

gene_means <- rowMeans( nrm_counts )
gene_vars <- rowVars_spm( nrm_counts )

plot( gene_means, gene_vars / gene_means,
      log = "xy", cex = .3, col = adjustcolor("black", alpha=.3), 
      xlab = "mean", ylab = "variance / mean" )

informative_genes <- names(which( 
  gene_vars / gene_means  >  3 * poisson_vmr ))


# variance stabilizing normalization
nrm_counts = as.matrix( t( sqrt( nrm_counts[ informative_genes,  ] ) ) )
