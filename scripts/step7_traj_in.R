suppressMessages(library("Seurat"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("slingshot"))
suppressMessages(library("SingleCellExperiment"))
suppressMessages(library("rgl"))
suppressMessages(library("rmarkdown"))
suppressMessages(library("gam"))
suppressMessages(library("pheatmap"))

# A. Parameters: folder configuration
dir.name = snakemake@params[["output_dir"]]
input_data = snakemake@input[["data"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis")

# B. Parameters: analysis configuration
selected_res = snakemake@params[["selected_res"]]
start.clus = snakemake@params[["start_clus"]]
end.clus = snakemake@params[["end_clus"]]
n_var_genes = snakemake@params[["n_var_genes"]]
n_plotted_genes = snakemake@params[["n_plotted_genes"]]


# C. Analysis
# 11. Trajectory inference analysis using slingshot
# 11.1 Seurat object to SingleCellExperiment object (required by slingshot) & cluster cluster_res 
seurat <- readRDS(input_data)
assay_type <- seurat@active.assay
cluster_res <- paste0(assay_type,"_snn_res.",selected_res)
if (!(cluster_res %in% colnames(seurat@meta.data))){
  stop("Specified resolution is not available.")
}
seurat.sim <- as.SingleCellExperiment(seurat)


# 11.2 Slingshot algorithm (dimensions = UMAP). There are 4 options depending of the start and end cluster in the following trajectory:
if(is.numeric(start.clus) == FALSE && is.numeric(end.clus) == FALSE){
  seurat.sim <- slingshot(seurat.sim, clusterLabels=cluster_res, reducedDim="UMAP")
  const = FALSE
} else if(is.numeric(start.clus) == TRUE && is.numeric(end.clus) == FALSE){
  seurat.sim <- slingshot(seurat.sim, clusterLabels=cluster_res, reducedDim="UMAP", start.clus = start.clus)
  const = TRUE
} else if(is.numeric(start.clus) == FALSE && is.numeric(end.clus) == TRUE){
  seurat.sim <- slingshot(seurat.sim, clusterLabels=cluster_res, reducedDim="UMAP", end.clus = end.clus)
  const = TRUE
} else if(is.numeric(start.clus) == TRUE && is.numeric(end.clus) == TRUE){
  seurat.sim <- slingshot(seurat.sim, clusterLabels=cluster_res, reducedDim="UMAP", start.clus = start.clus, end.clus=end.clus)
  const = TRUE
}

# 11.3 Plots

# 11.3.1 Set the number of colours from the palette and extend it.
n_col <- length(levels(seurat.sim@colData[, cluster_res]))
getPalette <- colorRampPalette(brewer.pal(9,'Set1'))

# 11.3.2 Set a good window size for the 3D plots.
r3dDefaults$windowRect <- c(0,50, 1024, 720) 

# 11.3.3 Curves 3D plot with legend --> HTML output.
plot3d(reducedDims(seurat.sim)$UMAP[,1:3], col = getPalette(n_col)[seurat.sim@colData[, cluster_res]])
plot3d(SlingshotDataSet(seurat.sim), lwd = 3, add = TRUE)
legend3d("topright", legend=paste0("Cluster - ", levels(seurat.sim@colData[, cluster_res])), pch=16, col=getPalette(n_col), inset=c(0.001))
writeWebGL(dir = paste0(dir.name, "/", folders[6]), filename = file.path(paste0(dir.name, "/", folders[6]), paste0("3D_curves_", selected_res, "_res.html")),  width = 1024)

# 11.3.4 Lineage 3D plot with legend --> HTML output.
plot3d(reducedDims(seurat.sim)$UMAP[,1:3], col = getPalette(n_col)[seurat.sim@colData[, cluster_res]]) 
plot3d(SlingshotDataSet(seurat.sim), lwd = 3, type = "lineages", add = TRUE, show.constraints = const)
legend3d("topright", legend=paste0("Cluster - ", levels(seurat.sim@colData[, cluster_res])), pch=16, col=getPalette(n_col), inset=c(0.001))
writeWebGL(dir = paste0(dir.name, "/", folders[6]), filename = file.path(paste0(dir.name, "/", folders[6]), paste0("3D_lineages_", selected_res, "_res.html")),  width = 1024)

# 11.3.5 Curves 2D plot with legend --> png output. 
png(paste0(dir.name, "/", folders[6], "/2D_curves_", selected_res, "_res.png"), width = 900, height = 900)
plot(reducedDims(seurat.sim)$UMAP, col = getPalette(n_col)[seurat.sim@colData[, cluster_res]],
     pch=16, asp = 1, main = "2D curves - Clusters Trajectories")
legend("topright", legend=paste0("Cluster - ", levels(seurat.sim@colData[, cluster_res])), pch=16, col=getPalette(n_col))
lines(SlingshotDataSet(seurat.sim), lwd=2, col = 'black')
dev.off()

# 11.3.6 Lineage 2D plot with legend --> png output.
png(paste0(dir.name, "/", folders[6], "/2D_lineages_", selected_res, "_res.png"), width = 900, height = 900)
plot(reducedDims(seurat.sim)$UMAP, col = getPalette(n_col)[seurat.sim@colData[, cluster_res]],
     pch=16, asp = 1, main = "2D lineages - Clusters Trajectories")
legend("topright", legend=paste0("Cluster - ", levels(seurat.sim@colData[, cluster_res])), pch=16, col=getPalette(n_col))
lines(SlingshotDataSet(seurat.sim), lwd=2, type = 'lineages', col = 'black', show.constraints = const)
dev.off()

# 11.4 Temporally expressed genes heatmap
# 11.4.1 Pseudotime is obtained
t <- seurat.sim$slingPseudotime_1

# 11.4.2 Get the n variable genes.
Y <- assays(seurat.sim)$logcounts
var_genes <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:n_var_genes]
Y <- Y[var_genes,]

# 11.4.3 Fitting a gam model
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  suppressWarnings({
    tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})

# 11.4.4 Get the top genes selected by p-val and plot the heatmap.
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:n_plotted_genes]
heatdata <- assays(seurat.sim)$logcounts[topgenes, order(t, na.last = NA)]
heatclus <- seurat.sim@colData[,cluster_res][order(t, na.last = NA)]
annotation <- data.frame("Cluster" = heatclus, row.names = colnames(heatdata))
ann_colors <- list("Cluster" = setNames(getPalette(n_col),
                                        levels(heatclus)))
# 11.4.5 Plot the genes.
pheatmap(heatdata, cluster_cols = FALSE, 
         color =  colorRampPalette(c("yellow", "red"))(100),
         annotation_col = annotation, annotation_colors = ann_colors,
         show_colnames = FALSE, filename = paste0(dir.name, "/", folders[6], "/Temporally_expressed_heatmaps_", selected_res, "_res.png"))

saveRDS(seurat.sim, file = paste0(dir.name, "/", folders[6], "/slingshot_sce.rds"))
