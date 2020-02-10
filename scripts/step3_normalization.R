suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))

# A. Parameters: folder configuration 
dir.name = snakemake@params[["output_dir"]]
input_data = snakemake@input[["data"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs")

# B. Parameters: analysis configuration 
normalization = snakemake@params[["normalization"]] # "sct" or "standard"
regress_out = snakemake@params[["regress_out"]] # true or false
vars_to_regress = c(snakemake@params[["vars_to_regress"]]) # check if null 
random_seed = snakemake@params[["random_seed"]]

# C. Analysis
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}
#Load seurat object 
seurat = readRDS(input_data)

# 5 Normalization
# 5.1 Normalize data depending of the method.
if(normalization == "standard"){
	seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
	seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2500)
	
	# Identify the 10 most highly variable genes
	top10 <- head(VariableFeatures(seurat), 15)
	p1 <- VariableFeaturePlot(seurat) + theme(legend.position="bottom") 
	LabelPoints(plot = p1, points = top10, repel = TRUE) + theme(legend.position="bottom") 
	ggsave(paste0(dir.name, "/",folders[1], "/6_variable_features.pdf"), plot = p1)
	
	# Scaling
	all.genes <- rownames(seurat)
	if(regress_out == TRUE){
		seurat <- ScaleData(seurat, features = all.genes, vars.to.regress = vars_to_regress)
	} else {
                seurat <- ScaleData(seurat, features = all.genes)			
	}
} else if (normalization == "SCT"){
	if(regress_out == TRUE){
                seurat <- SCTransform(seurat, vars.to.regress = vars_to_regress, verbose = FALSE)		
	} else {
                seurat <- SCTransform(seurat, verbose = FALSE)		
	}
} else {
	message("Normalization method not found.")
}

## 5.2 Run PCA
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), npcs = 50) # This result could all be saved in a table. 
# Visualizing PCA in Different Ways: elbow plot most variable genes 
p2 <- VizDimLoadings(seurat, dims = 1:2, reduction = "pca") + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[2], "/1_viz_dim_loadings.pdf"), plot = p2, scale = 1.5)#, height = height, width = height * aspect_ratio)
p3 <- DimPlot(seurat, reduction = "pca", pt.size = 0.5) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[2], "/2_dimplot.pdf"), plot = p3, scale = 1.5)
# 5.3. Determine the dimensionality of the dataset
seurat <- JackStraw(seurat, num.replicate = 100, dims = 50)
seurat <- ScoreJackStraw(seurat, dims = 1:50)
p4 <- ElbowPlot(seurat, ndims = 50) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[2], "/3_elbowplot.pdf"), plot = p4, scale = 1.5)
p5 <- JackStrawPlot(seurat, dims = 1:50) + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2, byrow=TRUE)) 
ggsave(paste0(dir.name, "/",folders[2], "/4_jackstrawplot.pdf"), plot = p5, scale = 2)

# Save RDS: we can use this object to generate all the rest of the data
saveRDS(seurat, file = paste0(dir.name, "/",folders[2], "/seurat_normalized-pcs.rds"))

