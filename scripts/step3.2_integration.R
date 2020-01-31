suppressMessages(library("Seurat"))
suppressMessages(library("ggplot2"))
suppressMessages(library("RColorBrewer"))

options(future.globals.maxSize = 2000*1024^2) 

# A. Parameters: folder configuration
dir.name = snakemake@params[["output_dir"]]
input_data = snakemake@input[["data"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis")
vars_to_regress = snakemake@params[["vars_to_regress"]]
norm_type = snakemake@params[["norm_type"]]
random_seed = snakemake@params[["random_seed"]]

# C. Analysis
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}

# 6 Integration.
# 6.1 Define a function to read each input files and perform the early steps of integration.
path_to_seurat_object <- function(x, vars_to_regress, norm_type) {
  experiment = head(tail(strsplit(x, "/")[[1]], n=3), n=1) #Assay name is stored to later use in integrated object metadata.
  seurat_obj <- readRDS(x)
  seurat_obj <- AddMetaData(object = seurat_obj,
                            metadata = experiment,
                            col.name = 'assay_name')
  if (norm_type == "SCT"){
    seurat_obj <- SCTransform(seurat_obj, vars.to.regress = vars_to_regress, return.only.var.genes = FALSE, verbose = FALSE)
  } else if (norm_type == "standard") {
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", 
                                       nfeatures = length(rownames(seurat_obj@assays$RNA@counts)), verbose = FALSE)
  }
  return(assign(paste0(experiment[[1]][6]),seurat_obj))
}

# 6.2 We apply the function obtaining a list with all the seurat objects.
seurat_object_list <- lapply(input_data, function(x) path_to_seurat_object(x, vars_to_regress, norm_type))

# 6.3 We find the minimum number of features between assays to get an aproximate feature number.
n_feat_calc <- function(seurat) {
  if (seurat@active.assay == "SCT"){
    return(length(rownames(seurat@assays$SCT@counts)))
  } else if (seurat@active.assay == "RNA"){
    return(length(rownames(seurat@assays$RNA@counts)))
  }
}
n_features <- min(unlist(lapply(seurat_object_list, function(x) n_feat_calc(x))))

# 6.4 Final integration steps and integration object creation.
if (norm_type == "SCT") {
  seurat.features <- SelectIntegrationFeatures(object.list = seurat_object_list, nfeatures = n_features, verbose = FALSE)
  seurat.list <- PrepSCTIntegration(object.list = seurat_object_list, anchor.features = seurat.features, 
                                    verbose = FALSE)
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT", 
                                           anchor.features = seurat.features, verbose = FALSE)
  seurat.integrated <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT", verbose = FALSE)
  
} else if (norm_type == "standard") {
  seurat.features <- SelectIntegrationFeatures(object.list = seurat_object_list, nfeatures = n_features, verbose = FALSE)
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat_object_list, dims = 1:100, anchor.features = seurat.features, verbose = FALSE)
  seurat.integrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:100, verbose = FALSE)
  seurat.integrated <- ScaleData(seurat.integrated, verbose = TRUE, vars.to.regress = vars_to_regress)  
}

# 6.5 PCA and Visualize Dimensional Reduction genes.
seurat.integrated <- RunPCA(seurat.integrated, ncps = 100, verbose = FALSE)
VizDimLoadings(seurat.integrated, dims = 1:2, reduction = "pca") + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[2], "/1_viz_dim_loadings.png"), scale = 1.5)

# 6.6 UMAP projection and integration visualization plot.
getPalette <- colorRampPalette(brewer.pal(9,'Set1'))
seurat.integrated <- RunUMAP(seurat.integrated, dims = 1:30, verbose = FALSE)
DimPlot(seurat.integrated, reduction = "umap", group.by = "assay_name", cols=getPalette(length(levels(as.factor(seurat.integrated$assay_name)))))
ggsave(paste0(dir.name, "/", folders[2], "/2_dimplot_UMAP.png"), plot = last_plot(), device = "png")

# 6.7 Principal component study using Elbow plot and Jack Straw Plot
seurat.integrated <- JackStraw(seurat.integrated, num.replicate = 100, dims = 30)
seurat.integrated <- ScoreJackStraw(seurat.integrated, dims = 1:30)
ElbowPlot(seurat.integrated, ndims = 30) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[2], "/3_elbowplot.png"), scale = 1.5)
JackStrawPlot(seurat.integrated, dims = 1:30) + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2, byrow=TRUE)) 
ggsave(paste0(dir.name, "/",folders[2], "/4_jackstrawplot.png"), scale = 2)

saveRDS(object = seurat.integrated, file = paste0(dir.name, "/", folders[2], "/seurat_normalized-pcs.rds"))
