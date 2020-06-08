log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("Seurat"))
suppressMessages(library("ggplot2"))
suppressMessages(library("RColorBrewer"))

# A. Parameters: folder configuration.
dir.name = snakemake@params[["output_dir"]]
input_data = snakemake@input[["data"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis")
vars_to_regress = snakemake@params[["vars_to_regress"]]
norm_type = snakemake@params[["norm_type"]]
random_seed = snakemake@params[["random_seed"]]
velocyto = snakemake@params[["velocyto"]]
outdir_config = snakemake@params[["outdir_config"]]
case = snakemake@params[["case"]]
ram = snakemake@res[["mem"]]

# C. Analysis.
options(future.globals.maxSize = ram*1024^2)

# Set sedd. 
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}

# Load cell cycle genes.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Set gene letter cases.
if (case == "lowercase") {
  s.genes <- str_to_lower(s.genes)
  g2m.genes <- str_to_lower(g2m.genes)
  message ("Set to lowercase.")
} else if (case == "titlecase") {
  s.genes <- str_to_title(s.genes)
  g2m.genes <- str_to_title(g2m.genes)
  message ("Set to titlecase.")
} else if (case == "uppercase") {
  message ("Set to uppercase.")
} else {
  message("Please choose a correc case option.")
  quit()
}


# 6. Integration.
# 6.1. Define a function to read each input files and perform the early steps of integration.
path_to_seurat_object <- function(x, vars_to_regress, norm_type, velocyto) {
  experiment = head(tail(strsplit(x, "/")[[1]], n=3), n=1) #Assay name is stored to later use in integrated object metadata.
  seurat_obj <- readRDS(x)
  seurat_obj <- AddMetaData(object = seurat_obj,
                            metadata = experiment,
                            col.name = 'assay_name')
  # If RNA velocity is going to be performed, we add the velocyto matrices in this step.
  if (velocyto){
    # Velocyto matrices path is infered. 
    velocyto_dir = paste0(outdir_config,"/star/", experiment,"/Solo.out/Velocyto/raw/")
    velo_names = c("spliced", "unspliced", "ambiguous")
    vel_matrices = list()
    # The matrices are read in 10X format.
    for (name in velo_names) {
      vel_matrices[[name]] <- Read10X(data.dir = paste0(velocyto_dir, name))
    }
    # The matrices are added as assays in the respective seurat object.
    for (name in velo_names) {
      vel_matrices[[name]] <- vel_matrices[[name]][, which(x = colnames(vel_matrices[[name]]) %in% colnames(seurat_obj))] 
      seurat_obj[[name]] <- CreateAssayObject(counts = vel_matrices[[name]])
    }
  }
  
  # The normalization is chosen.
  if (norm_type == "SCT"){
    seurat_obj <- SCTransform(seurat_obj, vars.to.regress = vars_to_regress, return.only.var.genes = FALSE, verbose = FALSE)
  } else if (norm_type == "standard") {
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", 
                                       nfeatures = length(rownames(seurat_obj@assays$RNA@counts)), verbose = FALSE)
  }
  return(assign(paste0(experiment[[1]][6]),seurat_obj))
}

# 6.2. The function is applied obtaining a list with all the seurat objects.
seurat_object_list <- lapply(input_data, function(x) path_to_seurat_object(x, vars_to_regress, norm_type, velocyto))

# 6.3. The minimum number of features between assays is obtained to get an aproximate feature number and keep as much genes as possible.
n_feat_calc <- function(seurat) {
  if (seurat@active.assay == "SCT"){
    return(length(rownames(seurat@assays$SCT@counts)))
  } else if (seurat@active.assay == "RNA"){
    return(length(rownames(seurat@assays$RNA@counts)))
  }
}
n_features <- min(unlist(lapply(seurat_object_list, function(x) n_feat_calc(x))))

# 6.4. Final integration steps and integration object creation.
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

# 6.5. PCA and Visualize Dimensional Reduction genes.
seurat.integrated <- RunPCA(seurat.integrated, ncps = 100, verbose = FALSE)
p1 <- VizDimLoadings(seurat.integrated, dims = 1:2, reduction = "pca") + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[2], "/1_viz_dim_loadings.pdf"), plot = p1, scale = 1.5)

# 6.6. PCA projection and integration visualization plot.
getPalette <- colorRampPalette(brewer.pal(9,'Set1'))
p2 <- DimPlot(seurat.integrated, reduction = "pca", group.by = "assay_name", cols=getPalette(length(levels(as.factor(seurat.integrated$assay_name)))))
ggsave(paste0(dir.name, "/", folders[2], "/2_dimplot_PCA.pdf"), plot = p2)

# 6.7. Principal component study using Elbow plot and Jack Straw Plot.
seurat.integrated <- JackStraw(seurat.integrated, num.replicate = 100, dims = 30)
seurat.integrated <- ScoreJackStraw(seurat.integrated, dims = 1:30)
p3 <- ElbowPlot(seurat.integrated, ndims = 30) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[2], "/3_elbowplot.pdf"), plot = p3, scale = 1.5)
p4 <- JackStrawPlot(seurat.integrated, dims = 1:30) + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2, byrow=TRUE)) 
ggsave(paste0(dir.name, "/",folders[2], "/4_jackstrawplot.pdf"), plot = p4, scale = 2)

# 6.8. Cell cycle scores and plots.
seurat.integrated <- CellCycleScoring(object = seurat.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
p5 <- FeaturePlot(object = seurat.integrated, features ="S.Score")
ggsave(paste0(dir.name, "/", folders[2], "/5_sscore_featureplot.pdf"), plot = p5, scale = 1.5)
p6 <- FeaturePlot(object = seurat.integrated, features ="G2M.Score")
ggsave(paste0(dir.name, "/", folders[2], "/6_g2mscore_featureplot.pdf"), plot = p6, scale = 1.5)
p7 <- DimPlot(seurat.integrated, reduction = "pca", pt.size = 0.5, label = TRUE, label.size = 5) + RotatedAxis() #+ theme(legend.position    ="bottom") 
ggsave(paste0(dir.name, "/", folders[2], "/7_no_umap_pca.pdf"), plot = p7, scale = 1.5)

# 6.9. Save expression matrix.
write.table(as.matrix(seurat.integrated@assays$integrated@scale.data), file = paste0(dir.name, "/", folders[2], "/normalized_expression_matrix.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# 6.10. Save Seurat object.
saveRDS(object = seurat.integrated, file = paste0(dir.name, "/", folders[2], "/seurat_normalized-pcs.rds"))
