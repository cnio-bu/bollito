library("Seurat")
library("dplyr")
library("data.table")
library("reticulate")
# library("clustree")
# py_install("umap-learn")  UMAP needs to be installed in order to run the RunUMAP function
library("ggplot2")

# A. Parameters: folder configuration 
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_celltypeid", "3_postprocessing", "4_degs", "5_gs")
cell_cycle_file = snakemake@params[["cc_file"]]
aspect_ratio <- 1.5
height <- 5
# B. Parameters: analysis configuration 
seed = snakemake@params[["seed"]]#randomly generate seed
pc = snakemake@params[["pc"]] # We should check the PCs using the Elbowplot
res = as.vector(snakemake@params[["res"]])

# C. Analysis
seurat <- readRDS(paste0(dir.name, "/", folders[2], "/seurat_normalized-pcs.rds"))
# 8. Post-processing
dir.create(paste0(dir.name, "/", folders[3]))
# 8.1. FindClusters
set.seed(seed)
seurat <- FindNeighbors(seurat, dims = 1:pc)
seurat <- FindClusters(seurat, resolution = res)
# 8.2 Clustree
#clustree(seurat, prefix = "RNA_snn_res.")
#ggsave(paste0(dir.name, "/", folders[3], "/1_Clustree.pdf"))
# 8.2. Run UMAP for all calculated resolutions
for(i in 1:length(which(grepl("RNA_snn_",colnames(seurat@meta.data))))){
  res = colnames(seurat@meta.data[which(grepl("RNA_snn_",colnames(seurat@meta.data)))][i])
  Idents(seurat) <- res
  seurat <- RunUMAP(seurat, dims = 1:50)
  DimPlot(seurat, reduction = "umap",  pt.size = 0.5, label = TRUE, label.size = 5) + theme_minimal() + theme(legend.position="bottom") 
  ggsave(paste0(dir.name, "/", folders[3], "/2_umap_",res,".pdf"))
}

FeaturePlot(seurat, 'nFeature_RNA', pt.size =  0.75) #+ theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/", folders[3], "/3_featureplot.pdf"))

# 8.3 Cell cycle
# Read in a list of cell cycle markers, from Tirosh et al, 2015.
# We can segregate this list into markers of G2/M phase and markers of S phase.
cc.genes <- readLines(cell_cycle_file)
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

seurat <- CellCycleScoring(object = seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
FeaturePlot(object = seurat, features ="S.Score") #+ theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/", folders[2], "/5_sscore_featureplot.pdf"))
FeaturePlot(object = seurat, features ="G2M.Score") #+ theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/", folders[2], "/6_g2mscore_featureplot.pdf"))

# 8.4. Visualize no - Umap
seurat.no.umap <- seurat
seurat.no.umap[["umap"]] <- NULL
DimPlot(seurat.no.umap, pt.size = 0.5, label = TRUE, label.size = 5) + RotatedAxis() #+ theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/", folders[3], "/4_no_umap.png"))


# Save RDS: we can use this object to generate all the rest of the data
saveRDS(seurat, file = paste0(dir.name, "/",folders[3], "/seurat_find-clusters.rds"))

