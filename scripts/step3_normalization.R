log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
suppressMessages(library("future"))

# A. Parameters: folder configuration. 
dir.name = snakemake@params[["output_dir"]]
input_data = snakemake@input[["seurat_obj"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs")

# B. Parameters: analysis configuration. 
normalization = snakemake@params[["norm_type"]] # "SCT" or "standard"
regress_out = snakemake@params[["regress_out"]] # true or false
vars_to_regress = c(snakemake@params[["vars_to_regress"]]) # check if null 
random_seed = snakemake@params[["random_seed"]]
regress_cell_cycle = snakemake@params[["regress_cell_cycle"]]
regress_merge_effect = snakemake@params[["regress_merge_effect"]]
case = snakemake@params[["case"]]
ram = snakemake@resources[["mem"]]
threads = snakemake@threads

# C. Analysis.
options(future.globals.maxSize = ram*1024^2)

# Set seed.
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}

#Set parallelization.
plan("multiprocess", workers = threads)

# Load cell cycle markers signature from Tirosh et al, 2015.
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

# Load seurat object.
seurat = readRDS(input_data)

# Regress merge variable input.
if (regress_merge_effect){
  merge_var = "assay_name"
} else {
  merge_var = NULL
}

# 5. Normalization.
# 5.1. Normalize data depending of the method.
if(normalization == "standard"){
  seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2500)

  # Identify the 10 most highly variable genes.
  top10 <- head(VariableFeatures(seurat), 10)
  p1 <- VariableFeaturePlot(seurat) + theme(legend.position="bottom") 
  LabelPoints(plot = p1, points = top10, repel = TRUE) + theme(legend.position="bottom") 
  ggsave(paste0(dir.name, "/",folders[1], "/6_variable_features.pdf"), plot = p1)
  
  # Scaling to perform PCA.
  # Merged object check.
  if (seurat@project.name == "merged") {
    seurat <- ScaleData(seurat, features = rownames(seurat), vars.to.regress = merge_var)
  } else { 
    seurat <- ScaleData(seurat, features = rownames(seurat))
  }
  # PCA previous to cell cycle scoring.
  seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), npcs = 50)

  # Cell cycle scores and plots.
  seurat <- CellCycleScoring(object = seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
  p4 <- FeaturePlot(object = seurat, features ="S.Score") + ggtitle("S phase score")
  ggsave(paste0(dir.name, "/", folders[2], "/4_sscore_featureplot.pdf"), plot = p4, scale = 1.5)
  p5 <- FeaturePlot(object = seurat, features ="G2M.Score") + ggtitle("G2/M phase score")
  ggsave(paste0(dir.name, "/", folders[2], "/5_g2mscore_featureplot.pdf"), plot = p5, scale = 1.5)
  p6 <- DimPlot(seurat, reduction = "pca", pt.size = 0.5, label = TRUE, label.size = 5) + RotatedAxis() #+ theme(legend.position    ="bottom") 
  ggsave(paste0(dir.name, "/", folders[2], "/6_cell_cycle_dimplot.pdf"), plot = p6, scale = 1.5)

  # Scaling.
  if(regress_out == TRUE){
    if (regress_cell_cycle) {
      if (seurat@project.name == "merged"){
        seurat <- ScaleData(seurat, features = rownames(seurat), vars.to.regress = c(vars_to_regress, "S.Score", "G2M.Score", merge_var))
      } else {
        seurat <- ScaleData(seurat, features = rownames(seurat), vars.to.regress = c(vars_to_regress, "S.Score", "G2M.Score"))
      }
      p7 <- DimPlot(seurat, group.by = "Phase", reduction = "pca", pt.size = 0.5, label = TRUE, label.size = 5) + RotatedAxis()
      ggsave(paste0(dir.name, "/", folders[2], "/6.1_cell_cycle_regressed_dimplot.pdf"), plot = p7, scale = 1.5)
    } else {
      if (seurat@project.name == "merged"){
        seurat <- ScaleData(seurat, features = rownames(seurat), vars.to.regress = c(vars_to_regress, merge_var))
      } else {
        seurat <- ScaleData(seurat, features = rownames(seurat), vars.to.regress = c(vars_to_regress))
      }      
    }
  }
} else if (normalization == "SCT"){
  if(regress_out == TRUE){
    if(seurat@project.name == "merged"){
      seurat <- SCTransform(seurat, vars.to.regress = c(vars_to_regress, merge_var), verbose = FALSE)
    } else {
      seurat <- SCTransform(seurat, vars.to.regress = vars_to_regress, verbose = FALSE)
    }
  } else {
    if(seurat@project.name == "merged"){
      seurat <- SCTransform(seurat, vars.to.regress = merge_var, verbose = FALSE)
    } else {
      seurat <- SCTransform(seurat, verbose = FALSE)
    }		
  }
  # PCA previous to cell cycle scoring.
  seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), npcs = 50) # This result could all be saved in a table.
  
  # Cell cycle scores and plots.
  seurat <- CellCycleScoring(object = seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  p4 <- FeaturePlot(object = seurat, features ="S.Score") + ggtitle("S phase score")
  ggsave(paste0(dir.name, "/", folders[2], "/4_sscore_featureplot.pdf"), plot = p4, scale = 1.5)
  p5 <- FeaturePlot(object = seurat, features ="G2M.Score") + ggtitle("G2/M phase score")
  ggsave(paste0(dir.name, "/", folders[2], "/5_g2mscore_featureplot.pdf"), plot = p5, scale = 1.5)
  p6 <- DimPlot(seurat, reduction = "pca", pt.size = 0.5, label = TRUE, label.size = 5) + RotatedAxis() #+ theme(legend.position    ="bottom") 
  ggsave(paste0(dir.name, "/", folders[2], "/6_cell_cycle_dimplot.pdf"), plot = p6, scale = 1.5)

  # If cell cycle regression is needed, a new SCT transformation is perform.
  if (regress_cell_cycle){
    if (regress_out) { 
      if (seurat@project.name == "merged"){
        seurat <- SCTransform(seurat, assay = "RNA", new.assay = "SCT", vars.to.regress = c(vars_to_regress, "S.Score", "G2M.Score", merge_var), verbose = FALSE)
      } else {
        seurat <- SCTransform(seurat, assay = "RNA", new.assay = "SCT", vars.to.regress = c(vars_to_regress, "S.Score", "G2M.Score"), verbose = FALSE)
      }
    } else {
      if (seurat@project.name == "merged"){
        seurat <- SCTransform(seurat, assay = "RNA", new.assay = "SCT", vars.to.regress = c("S.Score", "G2M.Score", merge_var), verbose = FALSE)
      } else {
        seurat <- SCTransform(seurat, assay = "RNA", new.assay = "SCT", vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
      }
    }
    p7 <- DimPlot(seurat, group.by = "Phase", reduction = "pca", pt.size = 0.5, label = TRUE, label.size = 5) + RotatedAxis()
    ggsave(paste0(dir.name, "/", folders[2], "/6.1_cell_cycle_regressed_dimplot.pdf"), plot = p7, scale = 1.5)
  }
} else {
	message("Normalization method not found.")
}

## 5.2. PCA metrics calculation.
Idents(seurat) <- seurat@project.name
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), npcs = 50) # This result could all be saved in a table. 
# Visualizing PCA in Different Ways: elbow plot most variable genes 
p2 <- VizDimLoadings(seurat, dims = 1:2, reduction = "pca") + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[2], "/1_viz_dim_loadings.pdf"), plot = p2, scale = 1.5)#, height = height, width = height * aspect_ratio)
p3 <- DimPlot(seurat, reduction = "pca", pt.size = 0.5) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[2], "/2_dimplot.pdf"), plot = p3, scale = 1.5)

# 5.3. Determine the dimensionality of the dataset
p4 <- ElbowPlot(seurat, ndims = 50) + theme(legend.position="bottom")
ggsave(paste0(dir.name, "/",folders[2], "/3_elbowplot.pdf"), plot = p4, scale = 1.5)

if (!(seurat@active.assay == "SCT")) {
  seurat <- JackStraw(seurat, num.replicate = 100, dims = 50)
  seurat <- ScoreJackStraw(seurat, dims = 1:50)
  p5 <- JackStrawPlot(seurat, dims = 1:50) + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2, byrow=TRUE)) + labs(y = "Empirical", x="Theoretical")
  ggsave(paste0(dir.name, "/",folders[2], "/4_jackstrawplot.pdf"), plot = p5, scale = 2)
}

if (seurat@project.name == "merged"){
  Idents(seurat) <- "assay_name"
  p6 <- DimPlot(seurat, reduction = "pca", pt.size = 0.5) + theme(legend.position="bottom")
  ggsave(paste0(dir.name, "/",folders[2], "/2_dimplot_merged.pdf"), plot = p6, scale = 1.5)
}

# 5.4. Save the expresion matrix.
if (normalization == "SCT") {
  write.table(as.matrix(seurat@assays$SCT@scale.data), file = paste0(dir.name, "/", folders[2], "/normalized_expression_matrix.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
}
if (normalization == "standard") {
  write.table(as.matrix(seurat@assays$RNA@scale.data), file = paste0(dir.name, "/", folders[2], "/normalized_expression_matrix.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
}

# 5.5. Save RDS: we can use this object to generate all the rest of the data.
saveRDS(seurat, file = paste0(dir.name, "/",folders[2], "/seurat_normalized-pcs.rds"))
