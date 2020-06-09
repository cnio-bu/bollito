log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))
suppressMessages(library("BiocParallel"))
suppressMessages(library("openxlsx"))

# A. Parameters: folder configuration. 
dir.name = snakemake@params[["output_dir"]]
input_data = snakemake@input[["seurat_obj"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis")

# B. Parameters: analysis configuration. 
selected_res = snakemake@params[["selected_res"]]
random_seed = snakemake@params[["random_seed"]]
test = snakemake@params[["test"]]

# C. Analysis.
# Set seed.
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}

# Load seurat object and set clustering resolution and assay type.
seurat <- readRDS(input_data)
assay_type <- seurat@active.assay
cluster_res <- paste0(assay_type, "_snn_res.", selected_res)
if (!(cluster_res %in% colnames(seurat@meta.data))){
  stop("Specified resolution is not available.")
}


#Set styles for xlsx files.
redStyle <- createStyle(fontColour = "#B60A1C", bgFill = "#FFF06A", textDecoration = c("BOLD"))
greenStyle <- createStyle(fontColour = "#309143", bgFill = "#FFF06A", textDecoration = c("BOLD"))


# 8. Differentially expressed genes between clusters. 
Idents(seurat) <- paste0(assay_type, "_snn_res.",selected_res)

# If the input seurat object is integrated we only perform the FindConservedMarkers.
if (seurat@active.assay == "integrated") {
  # 8.1. Table on differentially expressed genes - using basic filterings.
  for (i in 1:length(unique(Idents(seurat)))){
    clusterX.markers <- FindConservedMarkers(seurat, ident.1 = unique(Idents(seurat))[i], min.pct = 0.25, grouping.var = "assay_name", test.use = test) #min expressed
    write.table(clusterX.markers, file = paste0(dir.name, "/", folders[4], "/cluster", unique(Idents(seurat))[i],".markers.txt"), sep = "\t", col.names = NA, quote = FALSE)
  }
# If the input seurat object is not integrated.
} else {
  # 8.1. Table on differentially expressed genes - using basic filterings.
  for (i in 1:length(unique(Idents(seurat)))){
    clusterX.markers <- FindMarkers(seurat, ident.1 = unique(Idents(seurat))[i], min.pct = 0.25, test.use = test) #min expressed
    write.table(clusterX.markers, file = paste0(dir.name, "/", folders[4], "/cluster", unique(Idents(seurat))[i],".markers.txt"), sep = "\t", col.names = NA, quote = FALSE)
  }
  
  # 8.2. DE includying all genes - needed for a GSEA analysis. 
  for (i in 1:length(unique(Idents(seurat)))){
    clusterX.markers <- FindMarkers(seurat, ident.1 = unique(Idents(seurat))[i], min.pct = 0, logfc.threshold = 0, test.use = test) #min expressed
    #write.table(clusterX.markers, file = paste0(dir.name, "/", folders[4], "/cluster", unique(Idents(seurat))[i],".DE.txt"), sep = "\t", col.names = NA, quote = FALSE)
    wb <- createWorkbook()
    addWorksheet(wb, "DE analysis")
    writeData(wb, "DE analysis", clusterX.markers, rowNames = TRUE)
    conditionalFormatting(wb, "DE analysis", cols = 1:(ncol(clusterX.markers)+1),
                          rows = 2:(nrow(clusterX.markers) + 1), rule = "AND($C2<0, $F2<0.05)",
                          style = greenStyle)
    conditionalFormatting(wb, "DE analysis", cols = 1:(ncol(clusterX.markers)+1),
                          rows = 2:(nrow(clusterX.markers) + 1), rule = "AND($C2>0, $F2<0.05)",
                          style = redStyle)
    legend <- createComment(comment = c("Red means a positive LogFold\n\n", "Green means a negative LogFold"), style = c(redStyle, greenStyle))
    writeComment(wb, "DE analysis", col = 8, row = 2, comment = legend)
    saveWorkbook(wb, paste0(dir.name, "/", folders[4], "/cluster", unique(Idents(seurat))[i],".DE.xlsx"), overwrite = TRUE)    
    
    # 8.2.1. Create RNK file. 
    rnk = NULL
    rnk = as.matrix(clusterX.markers[,2])
    rownames(rnk)= toupper(row.names(clusterX.markers))
    write.table(rnk, file = paste0(dir.name, "/", folders[4], "/cluster", unique(Idents(seurat))[i],".rnk"), sep = "\t", col.names = FALSE, quote = FALSE)
  }
  
  # 8.3. Find TOP markers.
  seurat.markers <- FindAllMarkers(object = seurat, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25, test.use = test)
  groupedby.clusters.markers = seurat.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
  
  # 8.3.1. HeatMap top10.
  # setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
  p1 <- DoHeatmap(object = seurat, features = groupedby.clusters.markers$gene, cells = 1:length(colnames(seurat)), size = 8, angle = 45, 
	    group.bar = TRUE, draw.lines = F, raster = FALSE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + guides(color=FALSE) + theme(axis.text.y = element_text(size = 8)) + theme(legend.position="bottom") 
  ggsave(paste0(dir.name, "/", folders[4], "/1_heatmap_topmarkers.pdf"), plot = p1, scale = 3)
}

# 8.4. Save RDS: we can use this object to generate all the rest of the data.
saveRDS(seurat, file = paste0(dir.name, "/",folders[4], "/seurat_degs.rds"))
