log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
suppressMessages(library("Matrix"))

# A. Parameters: folder configuration 
data_dir = paste0(snakemake@params[["input_dir"]],"/Solo.out/Gene/filtered")
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis")

# B. Parameters: analysis configuration 
samples_path = snakemake@params[["samples_path"]]
min_cells_per_gene = snakemake@params[["min_cells_per_gene"]]
input_type = snakemake@params[["input_type"]]
units_path = snakemake@params[["units_path"]]
sample = snakemake@params[["sample"]]
random_seed = snakemake@params[["random_seed"]]

# C. Analysis
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}
# Read input and create the expression matrix object.
# If the input file is a fastq file (STARsolo input)
if (input_type == "fastq") {
  file.rename(paste0(data_dir,"/features.tsv"), paste0(data_dir,"/genes.tsv"))
  expression_matrix <- Read10X(data.dir = data_dir)
  rownames(expression_matrix) = stringr::str_to_title(rownames(expression_matrix))
#If the input file are matrices (directly read from units.tsv)
} else if (input_type == "matrix") { # units.tsv is loaded
  units <- read.csv(units_path, header = TRUE, sep = "\t", row.names = 1)
  if (units[sample,"unit"] == "10X") { #if input files are in 10X format (matrix, barcodes and genes files)
    expression_matrix <- readMM(toString(units[sample,"matrix"]))
    colnames(expression_matrix) <- read.table(toString(units[sample,"cell_names"]))[,1]
    row.names(expression_matrix) <- read.table(toString(units[sample,"genes"]))[,1]
  } else if (units[sample,"unit"] == "standard") { #if input file is a standard matris (genes as row names and cells as column names)
    expression_matrix = read.csv(toString(units[sample,"matrix"]), sep = "\t", header = TRUE, row.names = 1)
  } else {
    message("Please specify a correct unit input.")
  }
} else {
  message("Please specify a correct input type.")
}

# 1. Creating a seurat object 
seurat = CreateSeuratObject(expression_matrix, project = sample, min.features = 200, min.cells = min_cells_per_gene)

# 1.1 Add metadata  
samples_file = read.table(samples_path, sep = "\t", row.names = 1, header = TRUE)
for (i in 1:length(colnames(samples_file))) {
  seurat <- AddMetaData(seurat, samples_file[sample, i], col.name = colnames(samples_file)[i])
}

# 1.1.1 Add specific cell metadata.
if (input_type == "matrix") {
  if (file.exists(toString(units[sample,"metadata"]))){
    meta_file = read.table(toString(units[sample,"metadata"]), sep = "\t", row.names = 1, header = TRUE)
    meta_file <- subset(meta_file, row.names(meta_file) %in% row.names(seurat@meta.data))
    for (i in 1:length(colnames(meta_file))) {
        seurat <- AddMetaData(seurat, meta_file[,i], col.name = colnames(meta_file)[i])
    }
  }
} else {
  message("No metadata found.")
}

# 2. Preprocessing: Filter out low-quality cells
# 2.1. Mitochondrial genes - check levels of expression for mt genes 
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^Mt-")
# 2.2. Ribosomal genes - check levels of expression for rb genes 
seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^Rp[sl][[:digit:]]")
# 2.3. QC: violin plots
p1 <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.25) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/", folders[1], "/1_vlnplot_ngene_numi_pctmit_beforefilt.pdf"), plot = p1, scale = 1.5) 
p2 <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"), ncol = 3, pt.size = 0.25) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/", folders[1], "/2_vlnplot_ngene_numi_pctribo_beforefilt.pdf"), plot = p2, scale = 1.5)
# 2.4. QC: GenePlot
scatter1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.25)+ theme(legend.position="bottom") 
scatter2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.25) + theme(legend.position="bottom") 
p3 <- CombinePlots(plots = list(scatter1, scatter2))
ggsave(paste0(dir.name, "/", folders[1], "/3_geneplot_numi_vs_pctmit_ngene.pdf"), plot = p3, scale = 1.5)

# Save RDS: we can use this object to generate all the rest of the data
saveRDS(seurat, file = paste0(dir.name, "/" ,folders[1], "/seurat_pre-qc.rds"))

