library("Seurat")
library("dplyr")
library("data.table")
library("reticulate")
library("ggplot2")

# A. Parameters: folder configuration 
data_dir = paste0(snakemake@params[["input_dir"]],"/","Solo.out")
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_celltypeid", "3_postprocessing", "4_degs", "5_cellcycle")

# B. Parameters: analysis configuration 

# C. Analysis
seurat = readRDS(paste0(dir.name, "/", folders[1], "/seurat_post-qc.rds"))
# 5. Normalize data
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# 6. Find Variable Genes
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2500)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat), 15)
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
ggsave(paste0(dir.name, "/",folders[1], "/6_VariableFeaturesPlot.png"))


# 7. Start of Identifying Cell Types
dir.create(paste0(dir.name, "/", folders[2]))
# 7.1. Scale the data
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)
## 7.2. Run PCA
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), npcs = 100) # This result could all be saved in a table. 
# Visualizing PCA in Different Ways: elbow plot most variable genes 
VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
ggsave(paste0(dir.name, "/",folders[2], "/1_VizDimLoadings.png"))
DimPlot(seurat, reduction = "pca", pt.size = 1)
ggsave(paste0(dir.name, "/",folders[2], "/2_DimPlot.png"))
# 7.3. Determine the dimensionality of the dataset
seurat <- JackStraw(seurat, num.replicate = 100, dims = 100)
seurat <- ScoreJackStraw(seurat, dims = 1:100)
ElbowPlot(seurat, ndims = 100)
ggsave(paste0(dir.name, "/",folders[2], "/3_ElbowPlot.png"))
JackStrawPlot(seurat, dims = 1:100)
ggsave(paste0(dir.name, "/",folders[2], "/4_JackStrawPlot.png"))

# Save RDS: we can use this object to generate all the rest of the data
saveRDS(seurat, file = paste0(dir.name, "/",folders[2], "/seurat_normalized-pcs.rds"))

