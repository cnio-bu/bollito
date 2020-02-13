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

#Load cell cycle genes.
cc.genes <- c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Mlf1Ip","Hells","Rfc2","Rpa2","Nasp","Rad51Ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8Ap2","Usp1","Clspn","Pola1","Chaf1B","Brip1","E2F8","Hmgb2","Cdk1","Nusap1","Ube2C","Birc5","Tpx2","Top2A","Ndc80","Cks2","Nuf2","Cks1B","Mki67","Tmpo","Cenpf","Tacc3","Fam64A","Smc4","Ccnb2","Ckap2L","Ckap2","Aurkb","Bub1","Kif11","Anp32E","Tubb4B","Gtse1","Kif20B","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25C","Kif2C","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2E3","Gas2L3","Cbx5","Cenpa")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

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
p1 <- VizDimLoadings(seurat.integrated, dims = 1:2, reduction = "pca") + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[2], "/1_viz_dim_loadings.pdf"), plot = p1, scale = 1.5)

# 6.6 PCA projection and integration visualization plot.
getPalette <- colorRampPalette(brewer.pal(9,'Set1'))
p2 <- DimPlot(seurat.integrated, reduction = "pca", group.by = "assay_name", cols=getPalette(length(levels(as.factor(seurat.integrated$assay_name)))))
ggsave(paste0(dir.name, "/", folders[2], "/2_dimplot_PCA.pdf"), plot = p2)

# 6.7 Principal component study using Elbow plot and Jack Straw Plot
seurat.integrated <- JackStraw(seurat.integrated, num.replicate = 100, dims = 30)
seurat.integrated <- ScoreJackStraw(seurat.integrated, dims = 1:30)
p3 <- ElbowPlot(seurat.integrated, ndims = 30) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[2], "/3_elbowplot.pdf"), plot = p3, scale = 1.5)
p4 <- JackStrawPlot(seurat.integrated, dims = 1:30) + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2, byrow=TRUE)) 
ggsave(paste0(dir.name, "/",folders[2], "/4_jackstrawplot.pdf"), plot = p4, scale = 2)

# 6.8 Cell cycle scores and plots
seurat.integrated <- CellCycleScoring(object = seurat.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
p5 <- FeaturePlot(object = seurat.integrated, features ="S.Score")
ggsave(paste0(dir.name, "/", folders[2], "/5_sscore_featureplot.pdf"), plot = p5, scale = 1.5)
p6 <- FeaturePlot(object = seurat.integrated, features ="G2M.Score")
ggsave(paste0(dir.name, "/", folders[2], "/6_g2mscore_featureplot.pdf"), plot = p6, scale = 1.5)
p7 <- DimPlot(seurat.integrated, reduction = "pca", pt.size = 0.5, label = TRUE, label.size = 5) + RotatedAxis() #+ theme(legend.position    ="bottom") 
ggsave(paste0(dir.name, "/", folders[2], "/7_no_umap_pca.pdf"), plot = p7, scale = 1.5)

saveRDS(object = seurat.integrated, file = paste0(dir.name, "/", folders[2], "/seurat_normalized-pcs.rds"))
