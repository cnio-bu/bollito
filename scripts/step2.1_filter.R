log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

message("CONFIGURATION STEP")
# A. Parameters: 
# 1. Load libraries. 
suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))
message("1. Libraries were loaded.")

# 2. Folder configuration. 
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis")
message("2. Folder paths were set.")

# 3. Get variables from Snakemake. 
gene.filter = c(snakemake@params[["genes"]]) # Gene or list of genes.
filter.out = snakemake@params[["filter_out"]] # true or false
filter.threshold = snakemake@params[["threshold"]] # numeric
random_seed = snakemake@params[["random_seed"]]
ram = snakemake@resources[["mem"]]
message("3. Parameters were loaded.")

# 4. Analysis configuration. 
# RAM configuration.
options(future.globals.maxSize = ram*1024^2)
# Set seed.
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}
message(paste0("4. Seed was set at ", random_seed, "."))
message("Configuration finished.")
message("\n")

# B. Analysis.
message("PROCESSING STEP")
# Load Seurat object from previous step.
seurat = readRDS(paste0(dir.name, "/", folders[1], "/seurat_post-qc.rds"))
message("1. Seurat object was loaded.")

# 4. Filter cell based on gene expression values.
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
message("2. Cells were filter out.")

# Save RDS: we can use this object to generate all the rest of the data.
saveRDS(seurat, file = paste0(dir.name, "/",folders[1], "/seurat_post-qc-filtered.rds"))
message("3. Seurat object was saved.")


