suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))

# A. Parameters: folder configuration 
data_dir = paste0(snakemake@params[["input_dir"]],"/Solo.out/Gene/filtered")
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis")

# B. Parameters: analysis configuration 
project_name = snakemake@params[["project_name"]]
meta_path = snakemake@params[["meta_path"]]
min_cells_filter =  snakemake@params[["min_cells_filter"]]
# C. Analysis
# Read STARSolo output
file.rename(paste0(data_dir,"/features.tsv"), paste0(data_dir,"/genes.tsv"))
expression_matrix <- Read10X(data.dir = data_dir)
rownames(expression_matrix) = stringr::str_to_title(rownames(expression_matrix))

# Create Analysis folder
# 1. Creating a seurat object 
if (min_cells_filter) {
  seurat = CreateSeuratObject(expression_matrix, project = project_name, min.features = 200, min.cells = 3)
} else if (min_cells_filter == FALSE) {
  seurat = CreateSeuratObject(expression_matrix, project = project_name, min.features = 200)
} else {
  stop("Write True or False in config.yaml file please.")
}
# 1.1 Add metadata
metadata = read.csv(meta_path, sep = "\t", row.names = 1)
sample_name <- head(tail(unlist(strsplit(data_dir, "/")), n = 2), n= 1)
for (i in 1:length(colnames(metadata))) {
  seurat <- AddMetaData(seurat, metadata[sample_name, i], col.name = colnames(metadata)[i])
}

# 2. Preprocessing: Filter out low-quality cells
# 2.1. Mitochondrial genes - check levels of expression for mt genes 
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^Mt-")
# 2.2. Ribosomal genes - check levels of expression for rb genes 
seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^Rp[sl][[:digit:]]")
# 2.3. QC: violin plots
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.25) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/", folders[1], "/1_vlnplot_ngene_numi_pctmit_beforefilt.png"), scale = 1.5) 
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"), ncol = 3, pt.size = 0.25) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/", folders[1], "/2_vlnplot_ngene_numi_pctribo_beforefilt.png"), scale = 1.5)
# 2.4. QC: GenePlot
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.25)+ theme(legend.position="bottom") 
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.25) + theme(legend.position="bottom") 
CombinePlots(plots = list(plot1, plot2))
ggsave(paste0(dir.name, "/", folders[1], "/3_geneplot_numi_vs_pctmit_ngene.png"), scale = 1.5)

# Save RDS: we can use this object to generate all the rest of the data
saveRDS(seurat, file = paste0(dir.name, "/" ,folders[1], "/seurat_pre-qc.rds"))

