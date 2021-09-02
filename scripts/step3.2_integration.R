log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

message("CONFIGURATION STEP")
# A. Parameters: 
# 1. Load libraries. 
suppressMessages(library("Seurat"))
suppressMessages(library("ggplot2"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("stringr"))
suppressMessages(library("future"))
suppressMessages(library("stringr"))
message("1. Libraries were loaded.")

# 2. Folder configuration. 
dir.name = snakemake@params[["output_dir"]]
input_data = snakemake@input[["data"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis")
message("2. Folder paths were set.")

# 3. Get variables from Snakemake. 
regress_out = snakemake@params[["regress_out"]]
regress_cell_cycle = snakemake@params[["regress_cell_cycle"]]
vars_to_regress = snakemake@params[["vars_to_regress"]]
norm_type = snakemake@params[["norm_type"]]
vf = snakemake@params[["variable_features"]] 
random_seed = snakemake@params[["random_seed"]]
velocyto = snakemake@params[["velocyto"]]
outdir_config = snakemake@params[["outdir_config"]]
case = snakemake@params[["case"]]
ram = snakemake@resources[["mem"]]
threads = snakemake@threads
write_table = as.logical(snakemake@params[["write_table"]])
message("3. Parameters were loaded.")

# 4. Analysis configuration. 
# RAM configuration.
options(future.globals.maxSize = ram*1024^2)
#Set parallelization.
plan("multiprocess", workers = threads)
message(paste0("4. Threads were set at ", threads, "."))
# Set seed.
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}
message(paste0("5. Seed was set at ", random_seed, "."))

# 5. Load cell cycle markers signature from Tirosh et al, 2015.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

## 6. Change geneset signatures to specific gene case.
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
message(paste0("6. Case is set as ", case, "."))
message("Configuration finished.")
message("\n")


message("PROCESSING STEP")
# 6. Integration.
# 6.1. This function will read each one of the input files and will perform the early normalization steps for each sample.
path_to_seurat_object <- function(x, regress_out, vars_to_regress, regress_cell_cycle, norm_type, add_velocity) {
    experiment = head(tail(strsplit(x, "/")[[1]], n=3), n=1) #Assay name is stored to later use in integrated object metadata.
    seurat_obj <- readRDS(x)
    seurat_obj <- AddMetaData(object = seurat_obj, metadata = experiment, col.name = 'assay_name')
    # If RNA velocity is going to be performed, we add the velocyto matrices in this step.
    if (add_velocity){
   	# Velocyto matrices path is infered. 
	velocyto_dir = paste0(outdir_config,"/star/", experiment,"/Solo.out/Velocyto/raw/")
        velo_names = c("spliced", "unspliced", "ambiguous")
	vel_matrices = list()
	# The matrices are read in 10x format.
	for (name in velo_names) {
		vel_matrices[[name]] <- Read10X(data.dir = paste0(velocyto_dir, name))
	}
	# The matrices are added as assays in the respective seurat object.
	for (name in velo_names) {
		vel_matrices[[name]] <- vel_matrices[[name]][, which(x = colnames(vel_matrices[[name]]) %in% colnames(seurat_obj))] 
		seurat_obj[[name]] <- CreateAssayObject(counts = vel_matrices[[name]])
	}

    }
    # Save regression parameters variable
    regression_param <- c(regress_out, regress_cell_cycle)
    vars_param <- list(vars_to_regress = vars_to_regress, cell_cycle = c("S.Score", "G2M.Score"))

    # Normalization of the individual samples
    if(norm_type == "standard"){
        seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
        message(paste0("2. Seurat object was normalized using the ", norm_type, " approach"))
        seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2500)
        # Scaling to perform PCA.
        seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
        # PCA previous to cell cycle scoring.
        seurat_obj <- RunPCA(seurat_obj, features = if(vf) VariableFeatures(object = seurat_obj) else rownames(seurat_obj), npcs = 50)
        # Cell cycle scores and plots.
        seurat_obj <- CellCycleScoring(object = seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
        # Scaling.
        if(all(regression_param)){
            seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), vars.to.regress = c(vars_to_regress, "S.Score", "G2M.Score"))
        } else if (any(regression_param)){
            seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), vars.to.regress = unlist(vars_param[regression_param]))
        }
    } else if (norm_type == "SCT"){
        seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
	    message(paste0("2. Seurat object was normalized using the ", norm_type, " approach"))
        # PCA previous to cell cycle scoring.
        seurat_obj <- RunPCA(seurat_obj, features = if(vf) VariableFeatures(object = seurat_obj) else rownames(seurat_obj), npcs = 50)
        # Cell cycle scores and plots.
        seurat_obj <- CellCycleScoring(object = seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
        # If cell cycle regression is needed, a new SCT transformation is performed.
        if(all(regression_param)){
            seurat_obj <- SCTransform(seurat_obj, assay = "RNA", new.assay = "SCT", vars.to.regress = c(vars_to_regress, "S.Score", "G2M.Score"), verbose = FALSE)
        } else if (any(regression_param)){
            seurat_obj <- SCTransform(seurat_obj, assay = "RNA", new.assay = "SCT", vars.to.regress = unlist(vars_param[regression_param]), verbose = FALSE)
        }
    } else {
	    message("Normalization method not found.")
    }
    return(assign(paste0(experiment[[1]][6]),seurat_obj))
}

# 6.2. The function is applied obtaining a list with all the seurat objects.
seurat_object_list <- lapply(input_data, function(x) path_to_seurat_object(x, regress_out = regress_out, vars_to_regress = vars_to_regress, regress_cell_cycle = regress_cell_cycle, norm_type = norm_type, add_velocity = velocyto))
message("1. All Seurat objects were loaded.")

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
  seurat.list <- PrepSCTIntegration(object.list = seurat_object_list, anchor.features = seurat.features, verbose = FALSE)
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT", anchor.features = seurat.features, verbose = FALSE)
  seurat.integrated <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT", verbose = FALSE)
  message(paste0("2. All Seurat objects were integrated following the ", norm_type, " approach."))
} else if (norm_type == "standard") {
  seurat.features <- SelectIntegrationFeatures(object.list = seurat_object_list, nfeatures = n_features, verbose = FALSE)
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat_object_list, dims = 1:100, anchor.features = seurat.features, verbose = FALSE)
  seurat.integrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:100, verbose = FALSE)
  seurat.integrated <- ScaleData(seurat.integrated, verbose = TRUE)
  message(paste0("2. All Seurat objects were integrated following the ", norm_type, " approach."))
}

# 6.5. PCA and Visualize Dimensional Reduction genes.
seurat.integrated <- RunPCA(seurat.integrated, features = if(vf) VariableFeatures(object = seurat.integrated) else rownames(seurat.integrated), ncps = 100, verbose = FALSE)
p1 <- VizDimLoadings(seurat.integrated, dims = 1:2, reduction = "pca") + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[2], "/1_viz_dim_loadings.pdf"), plot = p1, scale = 1.5)

# 6.6. PCA projection and integration visualization plot.
getPalette <- colorRampPalette(brewer.pal(9,'Set1'))
p2 <- DimPlot(seurat.integrated, reduction = "pca", group.by = "assay_name", cols=getPalette(length(levels(as.factor(seurat.integrated$assay_name)))))
ggsave(paste0(dir.name, "/", folders[2], "/2_dimplot_PCA.pdf"), plot = p2)
message("3. PCA was done.")

# 6.7. Principal component study using Elbow plot.
p3 <- ElbowPlot(seurat.integrated, ndims = 50) + theme(legend.position="bottom")
ggsave(paste0(dir.name, "/",folders[2], "/3_elbowplot.pdf"), plot = p3, scale = 1.5)
message("4. Elbowplot was generated.")


# 6.8. Cell cycle scores and plots.
seurat.integrated <- CellCycleScoring(object = seurat.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
p5 <- FeaturePlot(object = seurat.integrated, features ="S.Score")
ggsave(paste0(dir.name, "/", folders[2], "/5_sscore_featureplot.pdf"), plot = p5, scale = 1.5)
p6 <- FeaturePlot(object = seurat.integrated, features ="G2M.Score")
ggsave(paste0(dir.name, "/", folders[2], "/6_g2mscore_featureplot.pdf"), plot = p6, scale = 1.5)
p7 <- DimPlot(seurat.integrated, reduction = "pca", pt.size = 0.5, label = TRUE, label.size = 5) + RotatedAxis() #+ theme(legend.position    ="bottom") 
ggsave(paste0(dir.name, "/", folders[2], "/7_no_umap_pca.pdf"), plot = p7, scale = 1.5)
message("6. Cell-cycle analysis plot was done.")

# 6.9. Save expression matrix.
if(write_table){
	write.table(as.matrix(seurat.integrated@assays$integrated@scale.data), file = paste0(dir.name, "/", folders[2], "/normalized_expression_matrix.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
	message("7. Integrated expression matrix was saved.")
} else {
	message("7. Integrated expression matrix was not saved, as specified in the configuration file.")
}

# 6.10. Save Seurat object.
saveRDS(object = seurat.integrated, file = paste0(dir.name, "/", folders[2], "/seurat_normalized-pcs.rds"))
message("8. Seurat object was saved.")
