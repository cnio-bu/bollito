suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))

# A. Parameters: folder configuration 
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_celltypeid", "3_postprocessing", "4_degs", "5_gs")
# B. Parameters: analysis configuration 
# C. Analysis
seurat = readRDS(paste0(dir.name, "/", folders[1], "/seurat_post-qc.rds"))
# 5. Normalize data
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# 6. Find Variable Genes
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2500)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat), 15)
plot1 <- VariableFeaturePlot(seurat) + theme(legend.position="bottom") 
LabelPoints(plot = plot1, points = top10, repel = TRUE) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[1], "/6_variable_features.pdf"))

# 7. Start of Identifying Cell Types
# 7.1. Scale the data
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)
## 7.2. Run PCA
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), npcs = 100) # This result could all be saved in a table. 
# Visualizing PCA in Different Ways: elbow plot most variable genes 
VizDimLoadings(seurat, dims = 1:2, reduction = "pca") + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[2], "/1_viz_dim_loadings.png"), scale = 1.5)#, height = height, width = height * aspect_ratio)
DimPlot(seurat, reduction = "pca", pt.size = 0.5) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[2], "/2_dimplot.png"), scale = 1.5)
# 7.3. Determine the dimensionality of the dataset
seurat <- JackStraw(seurat, num.replicate = 100, dims = 100)
seurat <- ScoreJackStraw(seurat, dims = 1:100)
ElbowPlot(seurat, ndims = 100) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[2], "/3_elbowplot.png"), scale = 1.5)
JackStrawPlot(seurat, dims = 1:100) + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2, byrow=TRUE)) 
ggsave(paste0(dir.name, "/",folders[2], "/4_jackstrawplot.png"), scale = 2)

# Save RDS: we can use this object to generate all the rest of the data
saveRDS(seurat, file = paste0(dir.name, "/",folders[2], "/seurat_normalized-pcs.rds"))

