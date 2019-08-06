suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))

# A. Parameters: folder configuration 
data_dir = paste0(snakemake@params[["input_dir"]],"/","Solo.out")
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_celltypeid", "3_postprocessing", "4_degs", "5_gs")

# B. Parameters: analysis configuration 
project_name = snakemake@params[["project_name"]]
# C. Analysis
# Read STARSolo output
expression_matrix <- Read10X(data.dir = data_dir)
rownames(expression_matrix) = stringr::str_to_title(rownames(expression_matrix))
# Create Analysis folder
# 1. Creating a seurat object 
seurat = CreateSeuratObject(expression_matrix, project = project_name, min.features = 200)

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

