NORMS = ["log", "SCT"]
PATH = "/icgc/dkfzlsdf/analysis/B210/Evelin/git-repo/"


#Create matrix.mtx, genes.tsv, barcodes.tsv for Seurat input	

#Make sure that all final output files get created
rule make_all:
	input:
		expand("{dataset}sce_RData/{norm}_decidua.RData", dataset=PATH, norm=NORMS)
	output:
		"out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done' > {output}"
	

#Create Seurat object from decidua data
rule decidua_normalization:
	input:
		"{dataset}decidua/matrix.mtx",
		"{dataset}decidua/barcodes.tsv",
		"{dataset}E-MTAB-6701_arrayexpress_10x_meta.txt"
	output:
		"{dataset}decidua/{norm}_seu.RData"
	shell:
		"""
		module load R/3.6.0 
		Rscript {wildcards.dataset}git-repo/decidua_normalization.R {wildcards.norm}
		"""

#Create SingleCellExperiment object
rule decidua_sce:
	input:
		"{dataset}decidua/{norm}_seu.RData"
	output: 
		"{dataset}sce_RData/{norm}_decidua.RData"
	shell:
		"""
		module load R/3.6.0 
		Rscript {wildcards.dataset}decidua_sce.R {wildcards.norm}
		"""		
