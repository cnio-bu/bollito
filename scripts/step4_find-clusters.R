suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("clustree"))
suppressMessages(library("ggplot2"))
suppressMessages(library("cluster"))
suppressMessages(library("writexl"))

# A. Parameters: folder configuration 
input_file = snakemake@input[["data"]]
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis")

# B. Parameters: analysis configuration 
pc = snakemake@params[["pc"]] # We should check the PCs using the Elbowplot
res = as.vector(snakemake@params[["res"]])
random_seed = snakemake@params[["random_seed"]]

# C. Analysis
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}
#Load seurat object
seurat <- readRDS(input_file)
assay_type <- seurat@active.assay

# 7. Post-processing
# 7.1 FindClusters using UMAP projection. We keep the significant PC obtained from PCA.
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:pc)
seurat <- FindClusters(seurat, resolution = res)
seurat <- RunUMAP(seurat,dims = 1:3, n.components = pc, verbose = FALSE)

# 7.2 Clustree
#clustree(seurat, prefix = "{assay_type}_snn_res.")
p1 <- clustree(seurat, prefix = paste0(assay_type,"_snn_res."))
ggsave(paste0(dir.name, "/", folders[3], "/1_Clustree.pdf"), plot = p1, scale = 1.5)

# 7.3 Clustering plots and silhouette parameters calculus.
# we create a empty lsit to store silhouette values.
silhouette_scores <- vector(mode = "list", length = length(res))

# loop for each resolution
for(i in 1:length(which(grepl(paste0(assay_type,"_snn_"),colnames(seurat@meta.data))))){
  full_res = colnames(seurat@meta.data[which(grepl(paste0(assay_type,"_snn_"),colnames(seurat@meta.data)))][i])
  Idents(seurat) <- full_res 
  p2 <- DimPlot(seurat, reduction = "umap", label = TRUE, label.size = 5) + theme_minimal() #+ theme(legend.position="bottom") 
  ggsave(paste0(dir.name, "/", folders[3], "/2_umap_",full_res,".pdf"), plot = p2, scale = 1.5)
 
  #Silhhouettes calculation
  dist.matrix <- dist(x = Embeddings(object = seurat[["pca"]])[, 1:pc])
  clusters <- slot(seurat, "meta.data")[,full_res]
  sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
  silhouette_scores[[i]] <- as.data.frame(summary(sil)[2])
  names(silhouette_scores[[i]]) <- full_res
}
#create a xlsx file to store the data.
write_xlsx(silhouette_scores, path = paste0(dir.name, "/", folders[3], "/3_silhouette_score.xlsx"),col_names = TRUE, format_headers = TRUE )

# 7.4 Feature plot
p3 <- FeaturePlot(seurat, 'nFeature_RNA', pt.size =  0.75) #+ theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/", folders[3], "/4_featureplot.pdf"), plot = p3, scale = 1.5)

# 7.5 Cell cycle
# Read in a list of cell cycle markers, from Tirosh et al, 2015.
# We can segregate this list into markers of G2/M phase and markers of S phase.
cc.genes <- c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Mlf1Ip","Hells","Rfc2","Rpa2","Nasp","Rad51Ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8Ap2","Usp1","Clspn","Pola1","Chaf1B","Brip1","E2F8","Hmgb2","Cdk1","Nusap1","Ube2C","Birc5","Tpx2","Top2A","Ndc80","Cks2","Nuf2","Cks1B","Mki67","Tmpo","Cenpf","Tacc3","Fam64A","Smc4","Ccnb2","Ckap2L","Ckap2","Aurkb","Bub1","Kif11","Anp32E","Tubb4B","Gtse1","Kif20B","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25C","Kif2C","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2E3","Gas2L3","Cbx5","Cenpa")
s.genes <- cc.genes[1:43] 
g2m.genes <- cc.genes[44:97]

seurat <- CellCycleScoring(object = seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
p4 <- FeaturePlot(object = seurat, features ="S.Score") + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/", folders[3], "/6_sscore_featureplot.pdf"), plot = p4, scale = 1.5)
p5 <- FeaturePlot(object = seurat, features ="G2M.Score") + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/", folders[3], "/7_g2mscore_featureplot.pdf"), plot = p5, scale = 1.5)

# 7.6 Visualize no - Umap
p6 <- DimPlot(seurat, reduction = "pca", pt.size = 0.5, label = TRUE, label.size = 5) + RotatedAxis() #+ theme(legend.position    ="bottom") 
ggsave(paste0(dir.name, "/", folders[3], "/5_no_umap_pca.pdf"), plot = p6, scale = 1.5)

# Save RDS: we can use this object to generate all the rest of the data
saveRDS(seurat, file = paste0(dir.name, "/",folders[3], "/seurat_find-clusters.rds"))

seurat
