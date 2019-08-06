configfile: "config.yaml"

#Extract variably expressed genes from decidua
rule extract_variable_genes:
	input:
		f1 = "{dataset}raw_data_10x.txt",
		f2 = "{dataset}E-MTAB-6701_arrayexpress_10x_meta.txt"
	output:
		f1 = "{dataset}seurat_object/decidua_5000_variably_expressed_counts.RData"
	resources:
		mem = 2000
	threads: 1
	shell:
		"""
		module load R/3.6.0
		Rscript --vanilla {wildcards.dataset}git-repo/decidua_pre.R 5000
		"""
		
#Extract variably expressed genes from decidua
rule extract_variable_genes:
	input:
		f1 = "{dataset}seurat_object/decidua_5000_variably_expressed_counts.RData",
		f3 = "{dataset}decidua_gene_list",
		f4 = "{dataset}seurat_object/joined.RData"
	output:
		f1 = "{dataset}seurat_object/decidua_sce.RData",
		f2 = "{dataset}seurat_object/menstrual_sce.RData"
	resources:
		mem = 2000
	threads: 1
	shell:
		"""
		module load R/3.6.0
		Rscript --vanilla {wildcards.dataset}git-repo/constructing_sce.R
		"""
		

#Make sure that all final output files get created
rule make_all:
	input:
		expand("{dataset}seurat_object/{type}_sce.RData", type=config["tissue"], dataset=config["path"])
	output:
		"out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done' > {output}"
		
