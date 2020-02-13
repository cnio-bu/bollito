log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))

# A. Parameters: folder configuration 
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis")


# B. Parameters: analysis configuration 
# project_name = "Test"
min_feat = snakemake@params[["min_feat"]]
max_feat = snakemake@params[["max_feat"]]
min_count = snakemake@params[["min_count"]]
max_count = snakemake@params[["max_count"]]
mit = snakemake@params[["mit"]]
ribo = snakemake@params[["ribo"]]
random_seed = snakemake@params[["random_seed"]]

# C. Analysis
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}
# Read RDS file from previous step
seurat = readRDS(paste0(dir.name, "/", folders[1], "/seurat_pre-qc.rds"))
# 3.1 We should apply the filterings once the QC plots (GenePlot and Violin plots) have been checked.
#seurat <- subset(seurat, subset = nFeature_RNA > min & nFeature_RNA < max & percent.mt < mit & percent.ribo < ribo )
# 3.1.1 Feature filter
if (!is.null(min_feat)) { 
  cells_seurat <- FetchData(object = seurat, vars = "nFeature_RNA") 
    if (is.null(max_feat)) {
      seurat <- seurat[, which(x = cells_seurat > min_feat)]                           
    } else {
      seurat <- seurat[, which(x = cells_seurat > min_feat & cells_seurat < max_feat)]
   }
}
# 3.1.2 Count filter
if (!is.null(min_count)) {
  cells_seurat <- FetchData(object = seurat, vars = "nCount_RNA")
  if (!is.null(max_feat)) {
    seurat <- seurat[, which(x = cells_seurat > min_count)] 
  } else {
    seurat <- seurat[, which(x = cells_seurat > min_count & cells_seurat < max_count)]
  }
}
# 3.1.3 Mitochondrial filter
if (!is.null(mit)) {
  mit_seurat <- FetchData(object = seurat, vars = "percent.mt")
  seurat <- seurat[, which(x = mit_seurat < mit)]
}
# 3.1.4 Ribosomal filter
if (!is.null(ribo)) {
ribo_seurat <- FetchData(object = seurat, vars = "percent.ribo")
seurat <- seurat[, which(x = ribo_seurat < ribo)]
}
# 3.2 QC: violin plots - After
p1 <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.25) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[1], "/4_vlnplot_ngene_numi_pctmit_afterfilt.pdf"), plot = p1, scale = 1.5)
p2 <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"), ncol = 3, pt.size = 0.25) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[1], "/5_vlnplot_ngene_numi_pctribo_afterfilt.pdf"), plot = p2, scale = 1.5)

# Save RDS: we can use this object to generate all the rest of the data
saveRDS(seurat, file = paste0(dir.name, "/",folders[1], "/seurat_post-qc.rds"))

