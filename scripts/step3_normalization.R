log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))

# A. Parameters: folder configuration 
dir.name = snakemake@params[["output_dir"]]
input_data = snakemake@input[["seurat_obj"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs")

# B. Parameters: analysis configuration 
normalization = snakemake@params[["normalization"]] # "sct" or "standard"
regress_out = snakemake@params[["regress_out"]] # true or false
vars_to_regress = c(snakemake@params[["vars_to_regress"]]) # check if null 
random_seed = snakemake@params[["random_seed"]]
regress_cell_cycle = snakemake@params[["regress_cell_cycle"]]

# C. Analysis
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}
#Load seurat object 
seurat = readRDS(input_data)

#Load cell cycle markers signature from Tirosh et al, 2015.
cc.genes <- c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Mlf1Ip","Hells","Rfc2","Rpa2","Nasp","Rad51Ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8Ap2","Usp1","Clspn","Pola1","Chaf1B","Brip1","E2F8","Hmgb2","Cdk1","Nusap1","Ube2C","Birc5","Tpx2","Top2A","Ndc80","Cks2","Nuf2","Cks1B","Mki67","Tmpo","Cenpf","Tacc3","Fam64A","Smc4","Ccnb2","Ckap2L","Ckap2","Aurkb","Bub1","Kif11","Anp32E","Tubb4B","Gtse1","Kif20B","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25C","Kif2C","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2E3","Gas2L3","Cbx5","Cenpa")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

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
  
  #Scaling to perform PCA.
  seurat <- ScaleData(seurat, features = rownames(seurat))

  # PCA previous to cell cycle scoring.
  seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), npcs = 50)

  # Cell cycle scores and plots.
  seurat <- CellCycleScoring(object = seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident 
= T)
  p4 <- FeaturePlot(object = seurat, features ="S.Score") + ggtitle("S phase score")
  ggsave(paste0(dir.name, "/", folders[2], "/4_sscore_featureplot.pdf"), plot = p4, scale = 1.5)
  p5 <- FeaturePlot(object = seurat, features ="G2M.Score") + ggtitle("G2/M phase score")
  ggsave(paste0(dir.name, "/", folders[2], "/5_g2mscore_featureplot.pdf"), plot = p5, scale = 1.5)
  p6 <- DimPlot(seurat, reduction = "pca", pt.size = 0.5, label = TRUE, label.size = 5) + RotatedAxis() #+ theme(legend.position    ="bottom") 
  ggsave(paste0(dir.name, "/", folders[2], "/6_cell_cycle_dimplot.pdf"), plot = p6, scale = 1.5)

  # Scaling
  if(regress_out == TRUE){
    if (regress_cell_cycle) {
      seurat <- ScaleData(seurat, features = rownames(seurat), vars.to.regress = c(vars_to_regress, "S.Score", "G2M.Score"))
    } else {
      seurat <- ScaleData(seurat, features = rownames(seurat), vars.to.regress = vars_to_regress)
    }
  }
} else if (normalization == "SCT"){
	if(regress_out == TRUE){
    seurat <- SCTransform(seurat, vars.to.regress = vars_to_regress, verbose = FALSE)		
	} else {
    seurat <- SCTransform(seurat, verbose = FALSE)		
	}
  #PCA previous to cell cycle scoring.
  seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), npcs = 50) # This result could all be saved in a table.
  
  #Cell cycle scores and plots.
  seurat <- CellCycleScoring(object = seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
  p4 <- FeaturePlot(object = seurat, features ="S.Score") + ggtitle("S phase score")
  ggsave(paste0(dir.name, "/", folders[2], "/4_sscore_featureplot.pdf"), plot = p4, scale = 1.5)
  p5 <- FeaturePlot(object = seurat, features ="G2M.Score") + ggtitle("G2/M phase score")
  ggsave(paste0(dir.name, "/", folders[2], "/5_g2mscore_featureplot.pdf"), plot = p5, scale = 1.5)
  p6 <- DimPlot(seurat, reduction = "pca", pt.size = 0.5, label = TRUE, label.size = 5) + RotatedAxis() #+ theme(legend.position    ="bottom") 
  ggsave(paste0(dir.name, "/", folders[2], "/6_cell_cycle_dimplot.pdf"), plot = p6, scale = 1.5)


  #If cell cycle regression is needed, a new SCT transformation is perform.
  if (regress_cell_cycle){
    seurat <- SCTransform(seurat, assay = "RNA", new.assay = "SCT", vars.to.regress = c(vars_to_regress, "S.Score", "G2M.Score"), verbose = FALSE)
  } 
} else {
	message("Normalization method not found.")
}

## 5.2 PCA metrics calculation.
Idents(seurat) <- seurat@project.name
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
p5 <- JackStrawPlot(seurat, dims = 1:50) + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2, byrow=TRUE)) + labs(y = "Empirical", x="Theoretical")
ggsave(paste0(dir.name, "/",folders[2], "/4_jackstrawplot.pdf"), plot = p5, scale = 2)

# Save RDS: we can use this object to generate all the rest of the data
saveRDS(seurat, file = paste0(dir.name, "/",folders[2], "/seurat_normalized-pcs.rds"))
