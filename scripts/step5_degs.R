library("Seurat")
library("dplyr")
library("data.table")
library("reticulate")
library("ggplot2")
library("BiocParallel")

# A. Parameters: folder configuration 
data_dir = paste0(snakemake@params[["input_dir"]],"/","Solo.out")
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_celltypeid", "3_postprocessing", "4_degs", "5_genesets")

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
scale_fill_gradientn(colors = c("blue", "white", "red")) + guides(color=FALSE) + theme(axis.text.y = element_text(size = 4))
ggsave(paste0(dir.name, "/", folders[4], "/1_Heatmap_TopMarkers_res0.1.pdf"))

# Number of identified genes by cluster
#n.clust <- length(unique(seurat@active.ident))
#markers <- lapply(seq_len(n.clust) - 1, function(cl) {
#	cl.markers <- FindMarkers(seurat, cl, logfc.threshold = 0, min.pct = 0.1, print.bar = FALSE)
#	cl.markers$cluster <- cl
#	cl.markers$gene <- rownames(cl.markers)
#	return(cl.markers)
#}) #, BPPARAM = bpparam)

#markers <- bind_rows(markers) %>%
#	select(gene, cluster, everything())

#markers.list <- lapply(0:(n.clust -1), function(x){
#	markers %>%
#		filter(cluster == x, pval<0.05) %>%
#		dplyr::arrange(-avg_logFC) %>%
#		select(Gene = gene, LogFC = avg_logFC, pVal = p_val)
#})
#names(markers.list) <- paste("Cluster", 0:(n.clust-1))

#marker.summary <- markers.list %>%
#	    map2_df(names(markers.list), ~ mutate(.x, Cluster = .y)) %>%
#	        mutate(IsUp = LogFC > 0) %>%
#		    group_by(Cluster) %>%
#		        summarise(Up = sum(IsUp), Down = sum(!IsUp)) %>%
#			    mutate(Down = -Down) %>%
#			        gather(key = "Direction", value = "Count", -Cluster) %>%
#				    mutate(Cluster = factor(Cluster, levels = names(markers.list)))

#ggplot(marker.summary,aes(x = fct_rev(Cluster), y = Count, fill = Direction)) +
#    geom_col() + geom_text(aes(y = Count + sign(Count) * max(abs(Count)) * 0.07, label = abs(Count)), size = 6, colour = "grey25") +
#    coord_flip() + scale_fill_manual(values = c("#377eb8", "#e41a1c")) + ggtitle("Number of identified genes") + 
#    theme(axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), legend.position = "bottom")

# Save RDS: we can use this object to generate all the rest of the data
saveRDS(seurat, file = paste0(dir.name, "/",folders[4], "/seurat_degs.rds"))
