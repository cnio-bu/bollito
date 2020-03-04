log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("clustree"))
suppressMessages(library("ggplot2"))
suppressMessages(library("cluster"))
suppressMessages(library("writexl"))

# A. Parameters: folder configuration 
input_file = snakemake@input[["seurat_obj"]]
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis")

# B. Parameters: analysis configuration 
pc = snakemake@params[["pc"]] # We should check the PCs using the Elbowplot
res = as.vector(snakemake@params[["res"]])
random_seed = snakemake@params[["random_seed"]]
k_neighbors = snakemake@params[["k_neighbors"]]

# C. Analysis
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}

#Load seurat object
seurat <- readRDS(input_file)
assay_type <- seurat@active.assay

# 7. Post-processing
# 7.1 FindClusters using UMAP projection. We keep the significant PC obtained from PCA.
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:pc, k.param = k_neighbors)
seurat <- FindClusters(seurat, resolution = res)
seurat <- RunUMAP(seurat,dims = 1:pc, n.components = 2, verbose = FALSE)

# 7.2 Clustree
#clustree(seurat, prefix = "{assay_type}_snn_res.")
p1 <- clustree(seurat, prefix = paste0(assay_type,"_snn_res."))
ggsave(paste0(dir.name, "/", folders[3], "/1_clustree.pdf"), plot = p1, scale = 1.5)

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
p3 <- FeaturePlot(seurat, 'nFeature_RNA', pt.size =  0.75) + labs(title = "Nº features") 
ggsave(paste0(dir.name, "/", folders[3], "/4_featureplot.pdf"), plot = p3, scale = 1.5)

# 7.5 Statistics table per cluster
# 7.5.1 Get the index of resolution columns.
resolutions = grep("snn_res", colnames(seurat@meta.data))
# 7.5.2 Set the highest number of clusters.
number_clusters = length(unique(seurat@meta.data$seurat_clusters))
# 7.5.3 Get the maximum cluster names vector.
cluster_names = paste0("Cluster ", seq(number_clusters)) 

# 7.5.4 Loop for each resolution and write table.
for(j in 1:length(resolutions)){
  # Set idents and levels.
  Idents(seurat) <- seurat@meta.data[,resolutions[j]]
  levels(Idents(seurat)) <- cluster_names[1:length(levels(Idents(seurat)))]
  table(Idents(seurat))
  # Proportion of cells per cluster
  prop.table(x = table(Idents(seurat)))
  # WhichCells(object = pbmc, ident = "0")
  cluster.averages <- AverageExpression(object = seurat)
  genes_per_cluster <- Matrix::colSums(cluster.averages$RNA>0)
  # Generate a joint table
  joint_stats = rbind(number_of_cells = table(Idents(seurat)), prop_per_cluster = prop.table(x = table(Idents(seurat))), genes_per_cluster = genes_per_cluster)
  colnames(joint_stats) = names(table(Idents(seurat)))
  write.table(t(joint_stats), file=paste0(dir.name, "/", folders[3], "/5_", colnames(seurat@meta.data)[resolutions[j]],"_stats.tsv"), sep="\t", col.names = NA, quote = FALSE)
}

# Save RDS: we can use this object to generate all the rest of the data
saveRDS(seurat, file = paste0(dir.name, "/",folders[3], "/seurat_find-clusters.rds"))
