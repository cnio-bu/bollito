library("Seurat")
library("dplyr")
library("data.table")
library("reticulate")

# A. Parameters: folder configuration 
data_dir = paste0(snakemake@params[["input_dir"]],"/","Solo.out")
dir.name = snakemake@params[["output_dir"]]
folders = c("1_Preprocessing", "2_CellTypeID", "3_Postprocessing", "4_DEGs", "5_Cell_cycle")

# B. Parameters: analysis configuration 
project_name = "Test"

# C. Analysis
# Read STARSolo output
expression_matrix <- Read10X(data.dir = data_dir)
# Create Analysis folder
dir.create(paste0(data_dir, dir.name))
# 1. Creating a seurat object 
seurat = CreateSeuratObject(expression_matrix, project = project_name, min.features = 200)

# 2. Preprocessing: Filter out low-quality cells
dir.create(paste0(dir.name, "/", folders[1]),recursive=TRUE)
# 2.1. Mitochondrial genes - check levels of expression for mt genes 
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
# 2.2. Ribosomal genes - check levels of expression for rb genes 
seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^Rp[sl][[:digit:]]")
# 2.3. QC: violin plots
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.5)
ggsave(paste0(dir.name, "/", folders[1], "/1_VlnPlot_distr_nGene_nUMI_percentMit_beforeFiltering.png"))
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"), ncol = 3, pt.size = 0.5)
ggsave(paste0(dir.name, "/", folders[1], "/2_VlnPlot_distr_nGene_nUMI_percentRibo_beforeFiltering.png"))
# 2.4. QC: GenePlot
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5)
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.5)
CombinePlots(plots = list(plot1, plot2))
ggsave(paste0(dir.name, "/", folders[1], "/3_GenePlot_nUMI_vs_percentMit_nGene.png"))


# Save RDS: we can use this object to generate all the rest of the data
saveRDS(seurat, file = paste0(paste0(dir.name, "/" ,folders[1], "/seurat_pre-QC.rds"))

