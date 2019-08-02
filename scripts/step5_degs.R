library("Seurat")
library("dplyr")
library("data.table")
library("reticulate")
library("ggplot2")
library("BiocParallel")

# A. Parameters: folder configuration 
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_celltypeid", "3_postprocessing", "4_degs", "5_gs")
# B. Parameters: analysis configuration 
selected_res = snakemake@params[["selected_res"]]

# C. Analysis
seurat <- readRDS(paste0(dir.name, "/", folders[3], "/seurat_find-clusters.rds"))
# 9 Differentially expressed genes between clusters. 
# dir.create(paste0(dir.name, "/", folders[4]))
# After checking out all results with the calculated resolutions, the rest of the analysis will be done using the specified one.
Idents(seurat) <- selected_res
# 9.1. Table on differentially expressed genes - using basic filterings
for (i in 1:length(unique(Idents(seurat)))){
  clusterX.markers <- FindMarkers(seurat, ident.1 = unique(Idents(seurat))[i], min.pct = 0.25) #min expressed
  print(x = head(x = clusterX.markers, n = 5))
  write.table(clusterX.markers, file=paste0(dir.name, "/", folders[4], "/cluster",unique(Idents(seurat))[i],".markers.txt"), sep="\t", col.names = NA)
}
# 9.2. DE includying all genes - needed for a GSEA analysis. 
for (i in 1:length(unique(Idents(seurat)))){
  clusterX.markers <- FindMarkers(seurat, ident.1 = unique(Idents(seurat))[i], min.pct = 0, logfc.threshold = 0) #min expressed
  print(x = head(x = clusterX.markers, n = 5))
  write.table(clusterX.markers, file=paste0(dir.name, "/", folders[4], "/cluster",unique(Idents(seurat))[i],".DE.txt"), sep="\t", col.names = NA)
  # Create RNK file 
  rnk = NULL
  rnk = as.matrix(clusterX.markers[,2])
  rownames(rnk)= toupper(row.names(clusterX.markers))
  write.table(rnk, file=paste0(dir.name, "/", folders[4], "/cluster",unique(Idents(seurat))[i],".rnk"), sep="\t", col.names = FALSE)
}
# 9.3. Find TOP markers
seurat.markers <- FindAllMarkers(object = seurat, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
groupedby.clusters.markers = seurat.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# HeatMap top10
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
DoHeatmap(object = seurat, features = groupedby.clusters.markers$gene, cells = 1:1000, size = 3, angle = 45, 
	  group.bar = TRUE, draw.lines = F, raster = FALSE) +
scale_fill_gradientn(colors = c("blue", "white", "red")) + guides(color=FALSE) + theme(axis.text.y = element_text(size = 4)) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/", folders[4], "/1_heatmap_topmarkers.pdf"), scale = 3)

# Save RDS: we can use this object to generate all the rest of the data
saveRDS(seurat, file = paste0(dir.name, "/",folders[4], "/seurat_degs.rds"))
