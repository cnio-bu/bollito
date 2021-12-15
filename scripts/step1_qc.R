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
suppressMessages(library("stringr"))
suppressMessages(library("Matrix"))
suppressMessages(library("patchwork"))
message("1. Libraries were loaded.")

# 2. Folder configuration. 
data_dir = paste0(snakemake@params[["input_dir"]],"/Solo.out/Gene/raw")
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis")
message("2. Folder paths were set.")

# 3. Get variables from Snakemake. 
samples_path = snakemake@params[["samples_path"]]
min_cells_per_gene = snakemake@params[["min_cells_per_gene"]]
input_type = snakemake@params[["input_type"]]
technology = snakemake@params[["technology"]]
units_path = snakemake@params[["units_path"]]
sample = snakemake@params[["sample"]]
random_seed = snakemake@params[["random_seed"]]
case =  snakemake@params[["case"]]
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

# 5. Set regexp for QC variables calculation.
if (case == "lowercase") {
  mito_grep <- "^mt-"
  ribo_grep <- "^rp[sl][[:digit:]]"
} else if (case == "titlecase") {
  mito_grep <- "^mt-"
  ribo_grep <- "^Rp[sl][[:digit:]]"
} else if (case == "uppercase") {
  mito_grep <- "^MT-"
  ribo_grep <- "^RP[SL][[:digit:]]"
} else {
  message("Please choose a correc case option.")
  quit()
}
message(paste0("5. Case is set as ", case, "."))
message("Configuration finished.")
message("\n")


message("PROCESSING STEP")
# 0. Read input and create the expression matrix object.
# If the input file is a fastq file (STARsolo input).
if (input_type == "fastq") {
  file.rename(paste0(data_dir,"/features.tsv"), paste0(data_dir,"/genes.tsv"))
  expression_matrix <- Read10X(data.dir = data_dir)
  message("1. Expression matrix was created from alignment files.")
# If the input file are matrices (directly read from units.tsv).
} else if (input_type == "matrix") { # units.tsv is loaded
  units <- read.csv(units_path, header = TRUE, sep = "\t", row.names = 1, comment.char = "#")
  # If the expression matrix is in 10x like format (matrix, cell barcodes and genes).
  if (technology == "10x") { 
    expression_matrix <- readMM(toString(units[sample,"matrix"]))
    colnames(expression_matrix) <- read.table(toString(units[sample,"cell_names"]))[,1]
    row.names(expression_matrix) <- read.table(toString(units[sample,"gene_names"]))[,1]
    message("1. Expression matrix was created from 10x/CellRanger-like files.")
  # If the expression matrix is in TSV format ((genes as row names and cells as column names).
  } else if (technology == "standard") { 
    expression_matrix = read.csv(toString(units[sample,"matrix"]), sep = "\t", header = TRUE, row.names = 1)
    message("1. Expression matrix was created from TSV matrix.")

  } else {
    message("Please specify a correct unit input.")
  }
} else {
  message("Please specify a correct input type.")
}

# 1. Creating a seurat object. 
expression_matrix <- expression_matrix[, colSums(expression_matrix != 0) > 100] # Take into account cells with more than 100 counts, since CreateSeuratObject function breaks.
seurat = CreateSeuratObject(expression_matrix, project = sample, min.features = 200, min.cells = min_cells_per_gene)
message("2. Seurat object was created.")

# 1.1 Add metadata from samples.tsv file.  
samples_file = read.table(samples_path, sep = "\t", row.names = 1, header = TRUE)
if (dim(samples_file)[2] != 0){
  for (i in 1:length(colnames(samples_file))) {
    seurat <- AddMetaData(seurat, samples_file[sample, i], col.name = colnames(samples_file)[i])
  }
}
message("3. Metadata from sample.tsv was added to Seurat object.")

# 1.1.1 Add specific cell metadata from metadata.tsv file.
if (input_type == "matrix") {
  if (file.exists(toString(units[sample,"metadata"]))){
    meta_file = read.table(toString(units[sample,"metadata"]), sep = "\t", row.names = 1, header = TRUE)
    meta_file <- subset(meta_file, row.names(meta_file) %in% row.names(seurat@meta.data))
    for (i in 1:length(colnames(meta_file))) {
        seurat <- AddMetaData(seurat, meta_file[,i], col.name = colnames(meta_file)[i])
    }
    message("3.1 Metadata from metadata.tsv was added to Seurat object.")
  } else {
    message("Metadata.tsv was not found.")
  }
}

#1.2 Set idents to avoid new idents based on shared CB names. 
Idents(seurat) <- rep(sample, length(colnames(seurat$RNA@data)))

# 2. Preprocessing: Get QC-related values.
# 2.1. Mitochondrial genes - check levels of expression for mt genes. 
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = mito_grep)
# 2.2. Ribosomal genes - check levels of expression for rb genes.
seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = ribo_grep)
message("4. Ribosomal and mitochondrial percentages were calculated.")

# 2.3. QC: violin plots.
p1 <- VlnPlot(seurat, features = c("nFeature_RNA"), pt.size = 0.25, cols = "#9CCCD0") + ggtitle("Nº features") + theme(legend.position="bottom") 
p2 <- VlnPlot(seurat, features = c("nCount_RNA"), pt.size = 0.25, cols = "#8ADD56")  + ggtitle("Nº counts") + theme(legend.position="bottom")
p3 <- VlnPlot(seurat, features = c("percent.mt"), pt.size = 0.25, cols = "#F07800") + ggtitle("Mitochondrial %") + theme(legend.position="bottom")
p4 <- VlnPlot(seurat, features = c("percent.ribo"), pt.size = 0.25, cols = "#E44631") + ggtitle("Ribosomal %") + theme(legend.position="bottom")
p_comp <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4)
ggsave(paste0(dir.name, "/", folders[1], "/1_vlnplot_QC_variables_prefilt.pdf"), plot = p_comp, scale = 1.2, width = 10, height = 8)
message("5. Combined violin plot was generated.")

# 2.4. QC: GenePlot.
scatter1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.25)+ theme(legend.position="bottom") + labs(title = "Mitochondrial % vs Nº counts", x = "Nº counts", y = "Mitochondrial %")
scatter2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.25) + theme(legend.position="bottom") + labs(title = "Nº features vs Nº counts", x = "Nº counts", y = "Nº features")
p_comb2 <- scatter1 + scatter2 + plot_layout(ncol = 2)
ggsave(paste0(dir.name, "/", folders[1], "/2_geneplot_numi_vs_pctmit_ngene.pdf"), plot = p_comb2, scale = 1.5)
message("6. Scaterplots were generated.")


# 2.5. Save RDS: we can use this object to generate all the rest of the data.
saveRDS(seurat, file = paste0(dir.name, "/" ,folders[1], "/seurat_pre-qc.rds"))
message("7. Seurat object was saved.")
