log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))

# A. Parameters: folder configuration.
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis")

# B. Parameters: analysis configuration. 
gene.filter = c(snakemake@params[["gene"]]) # Gene or list of genes.
filter.out = snakemake@params[["filter_out"]] # true or false
filter.threshold = snakemake@params[["threshold"]] # numeric
random_seed = snakemake@params[["random_seed"]]

# C. Analysis.
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}
# Read RDS file from previous step.
seurat = readRDS(paste0(dir.name, "/", folders[1], "/seurat_post-qc.rds"))

# 4.1 If there are negative markers availale: filter out cells based on gene expression. In this specific case, we are filtering out all cells expressing: Epcam, Pecam1, Krt19 and Ptprc. CHECK THIS

if(length(gene.filter) > 0){
	if (filter.out == TRUE){
		for(i in 1:length(gene.filter)){
			sub_seurat <- FetchData(object = seurat, vars = gene.filter[i])
			seurat <- seurat[, which(x = sub_seurat > filter.threshold)]
		}			
	} else {
		for(i in 1:length(gene.filter)){
			sub_seurat <- FetchData(object = seurat, vars = gene.filter[i])
			seurat <- seurat[, which(x = sub_seurat <= filter.threshold)]
		}	
	}
} else {
	next()
}

# Save RDS: we can use this object to generate all the rest of the data.
saveRDS(seurat, file = paste0(dir.name, "/",folders[1], "/seurat_post-qc-filtered.rds"))

