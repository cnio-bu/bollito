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
suppressMessages(library("patchwork"))
message("1. Libraries were loaded.")

# 2. Folder configuration. 
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis")
message("2. Folder paths were set.")

# 3. Get variables from Snakemake. 
min_feat = snakemake@params[["min_feat"]]
max_feat = snakemake@params[["max_feat"]]
min_count = snakemake@params[["min_count"]]
max_count = snakemake@params[["max_count"]]
mit = snakemake@params[["mit"]]
ribo = snakemake@params[["ribo"]]
random_seed = snakemake@params[["random_seed"]]
ram = snakemake@resources[["mem"]]
write_table = as.logical(snakemake@params[["write_table"]])
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
seurat = readRDS(paste0(dir.name, "/", folders[1], "/seurat_pre-qc.rds"))
message("1. Seurat object was loaded.")

#3. Preprocessing: Filter out low-quality cells.
# 3.1 Pre-filtering stats calculus.
stats_pre <- c(length(colnames(seurat)), median(seurat@meta.data[["nCount_RNA"]]), median(seurat@meta.data[["nFeature_RNA"]]), median(seurat@meta.data[["percent.mt"]]), median(seurat@meta.data[["percent.ribo"]]))
message("2. Pre-filtering statistics were obtained.")


# 3.2 We should apply the filterings once the QC plots (GenePlot and Violin plots) have been checked.
# 3.2.1 Feature filter.
min_feat <- ifelse(is.null(min_feat), 0, min_feat)
max_feat <- ifelse(is.null(max_feat), 0, max_feat)
if (min_feat != 0 | max_feat != 0) { 
  cells_seurat <- FetchData(object = seurat, vars = "nFeature_RNA")
  seurat <- seurat[, which(x = cells_seurat > min_feat & cells_seurat < max_feat)]
}

# 3.2.2 Count filter.
min_count <- ifelse(is.null(min_count), 0, min_count)
max_count <- ifelse(is.null(max_count), 0, max_count)
if (min_count != 0 | max_count != 0) { 
  cells_seurat <- FetchData(object = seurat, vars = "nCount_RNA")
  seurat <- seurat[, which(x = cells_seurat > min_count & cells_seurat < max_count)]
}


# 3.2.3 Mitochondrial filter.
if (!is.null(mit)) {
  mit_seurat <- FetchData(object = seurat, vars = "percent.mt")
  seurat <- seurat[, which(x = mit_seurat < mit)]
}

# 3.2.4 Ribosomal filter.
if (!is.null(ribo)) {
  ribo_seurat <- FetchData(object = seurat, vars = "percent.ribo")
  seurat <- seurat[, which(x = ribo_seurat < ribo)]
}
message("3. Cell filters were applied.")

# 3.3 QC: violin plots - After filter.
p1 <- VlnPlot(seurat, features = c("nFeature_RNA"), pt.size = 0.25, cols = "#9CCCD0") + ggtitle("Nº features") + theme(legend.position="bottom") 
p2 <- VlnPlot(seurat, features = c("nCount_RNA"), pt.size = 0.25, cols = "#8ADD56")  + ggtitle("Nº counts") + theme(legend.position="bottom")
p3 <- VlnPlot(seurat, features = c("percent.mt"), pt.size = 0.25, cols = "#F07800") + ggtitle("Mitochondrial %") + theme(legend.position="bottom")
p4 <- VlnPlot(seurat, features = c("percent.ribo"), pt.size = 0.25, cols = "#E44631") + ggtitle("Ribosomal %") + theme(legend.position="bottom")
p_comp <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4)
ggsave(paste0(dir.name, "/", folders[1], "/3_vlnplot_QC_variables_postfilt.pdf"), plot = p_comp, scale = 1.2, width = 10, height = 8)
message("4. Combined violin plot post-filtering was generated.")

# 3.4 Post-filter stats calculus.
stats_post <- c(length(colnames(seurat)), median(seurat@meta.data[["nCount_RNA"]]), median(seurat@meta.data[["nFeature_RNA"]]), median(seurat@meta.data[["percent.mt"]]), median(seurat@meta.data[["percent.ribo"]]))
message("5. Post-filtering statistics were obtained.")

# 3.5 Statistics table.
filtering_df <- data.frame("Number of cells" = c(stats_pre[1],stats_post[1]), "Count median" = c(stats_pre[2],stats_post[2]),"Expressed genes median" = c(stats_pre[3],stats_post[3]), "Mitochondrial percentage median" = c(stats_pre[4],stats_post[4]), "Ribosomal percentage median" = c(stats_pre[5],stats_post[5]))
row.names(filtering_df) <- c("Pre-QC", "Post-QC")
write.table(filtering_df, file = paste0(dir.name, "/", folders[1], "/4_pre_vs_post_stats.tsv"), sep = "\t", col.names = NA, quote = FALSE)
message("6. Statistics table was saved.")

# 3.6 Save expression matrix.
if(write_table){
	write.table(as.matrix(seurat@assays$RNA@counts), file = paste0(dir.name, "/", folders[1], "/expression_matrix.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
	message("7. Post-qc expression matrix was saved.")
} else {
	message("7. Post-qc expression matrix was not saved, as specified in the configuration file.")
}

# 3.7 Save RDS: we can use this object to generate all the rest of the data.
saveRDS(seurat, file = paste0(dir.name, "/",folders[1], "/seurat_post-qc.rds"))
message("8. Seurat object was saved.")
