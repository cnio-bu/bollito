log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))
suppressMessages(library("qusage"))
suppressMessages(library("clustree"))
suppressMessages(library("patchwork"))

# A. Parameters: folder configuration 
dir.name = snakemake@params[["output_dir"]]
input_data = snakemake@input[["seurat_obj"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis")

# B. Parameters: analysis configuration 
geneset_collection = snakemake@params[["gs_collection"]]
random_seed = snakemake@params[["random_seed"]]
resolutions = snakemake@params[["resolutions"]]
geneset_percentage <- snakemake@params[["geneset_percentage"]]

# C. Analysis
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}
# Load seurat object
seurat <- readRDS(input_data)
assay_type <- seurat@active.assay
dir.create(paste0(dir.name, "/", folders[5]))

# 9. GS scoring
# 9.1 Geneset loading, filtering and scoring.
genesets <- read.gmt(geneset_collection) #should be a tab file, each row = pathway.
genesets <- genesets[unlist(lapply(genesets, function(x) ((length(which(x%in% rownames(seurat)))/ length(x))> geneset_percentage)))]
seurat <- AddModuleScore(object = seurat, features= genesets, name = names(genesets))


# 9.2 We create a vector containing the different combinations for each resolutions and each cluster.  
res_clus_comb <- vector()  #Empty vector.
for (res in resolutions){
  full_res = paste0(assay_type, "_snn_res.", res)
  for (cluster in levels(seurat@meta.data[[full_res]])) {
    x <- paste0(full_res,"C",cluster) 
    res_clus_comb <- append(res_clus_comb, x) #Adding combinations to the vector.
  }
}

# 9.3 We create a MxN matrix for the p-values, where M is the previous calculated combinations and N is the different gene sets used.
mtx_pval <- matrix(nrow=length(res_clus_comb), ncol=length(names(genesets)))


# 9.4 We fill the matrix with the p-value calculated at cluster level for each resolution and gene set.
for (col in 1:length(genesets)) {
  row = 1
  for (res in resolutions){
    full_res = paste0(assay_type, "_snn_res.", res)
    for (cluster in levels(seurat@meta.data[[full_res]])) {
      # The resolution is set as Idents
      Idents(seurat) <- seurat[[full_res]]
      
      #Counts per gene calculus.
      seurat_target_clust <- subset(seurat, idents = cluster)
      counts_in_gene_set_target<- rowMeans(as.matrix(seurat_target_clust@assays[[assay_type]]@counts))[which(names(rowSums(as.matrix(seurat_target_clust@assays[[assay_type]]@counts))) %in% genesets[[col]])]
      #expr_genes_target <- length(counts_in_gene_set_target[counts_in_gene_set_target > 1])
      #all_expr_genes_target <- length(rowMeans(as.matrix(seurat_target_clust@assays[[assay_type]]@counts))[rowSums(as.matrix(seurat_target_clust@assays[[assay_type]]@counts)) > 1])
      
      seurat_offtarget_clust <- subset(seurat, idents = cluster, invert = TRUE)
      counts_in_gene_set_offtarget<- rowMeans(as.matrix(seurat_offtarget_clust@assays[[assay_type]]@counts))[which(names(rowSums(as.matrix(seurat_offtarget_clust@assays[[assay_type]]@counts))) %in% genesets[[col]])]
      #expr_genes_offtarget <- length(counts_in_gene_set_offtarget[counts_in_gene_set_offtarget > 1])
      #all_expr_genes_offtarget <- length(rowMeans(as.matrix(seurat_offtarget_clust@assays[[assay_type]]@counts))[rowSums(as.matrix(seurat_offtarget_clust@assays[[assay_type]]@counts)) > 1])
      
      #Max count value used to generate the factor. 
      max_val <- max(c(counts_in_gene_set_target, counts_in_gene_set_offtarget))
      #gene_matrix <- matrix(data = c(expr_genes_target, expr_genes_offtarget, all_expr_genes_target, all_expr_genes_offtarget), nrow = 2)    
      
      #The factor is calculated and applied for each vector.
      target_factor <- length(counts_in_gene_set_target[counts_in_gene_set_target > 0.05*max_val])/length(counts_in_gene_set_target)
      offtarget_factor <- length(counts_in_gene_set_offtarget[counts_in_gene_set_offtarget > 0.05*max_val])/length(counts_in_gene_set_offtarget)
      
      counts_in_gene_set_target_factor <- counts_in_gene_set_target*target_factor
      counts_in_gene_set_offtarget_factor <- counts_in_gene_set_offtarget*offtarget_factor
      
      # Wilcoxon test: paired, using a vectors of gene means for each Seurat subset.
      p_value <- wilcox.test(counts_in_gene_set_target_factor, counts_in_gene_set_offtarget_factor, paired = TRUE, alternative = "greater")$p.value
      print(p_value)
      mtx_pval[row, col] <- p_value
      row = row + 1
    }
  }
}

# 9.5 The p-values are corrected oer geneset.
mtx_corr_pval <- mtx_pval
for (col in 1:length(genesets)) {
  p_val_vec_corr <- p.adjust(mtx_corr_pval[,col], "fdr") 
  mtx_corr_pval[,col] <- p_val_vec_corr
}

# 9.6 We create a dataframe from a matrix
df_pval <- data.frame(mtx_corr_pval)
row.names(df_pval) <- res_clus_comb
colnames(df_pval) <- names(genesets)
df_pval <- df_pval[ order(row.names(df_pval)), ]

# 9.7 Save p-value table
write.csv(df_pval, file = paste0(dir.name, "/",folders[5], "/pval_table.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# 9.8 We loop for each geneset generating the plots
for (i in 1:length(genesets)){
  gs_name = paste0(names(genesets[i]), i)
  module_name = colnames(seurat@meta.data)[grep(names(genesets)[i], colnames(seurat@meta.data))]
  
  #Feature plot
  p1 <- FeaturePlot(object = seurat, features = module_name) #+ theme(legend.position="bottom") 
  ggsave(paste0(dir.name, "/", folders[5], "/", names(genesets)[i], "_featureplot.pdf"), plot = p1, scale = 1.5)
  
  #Mean expression plot
  p2 <- clustree(seurat, paste0(assay_type,"_snn_res."), node_colour = gs_name , node_colour_aggr = "mean") + ggtitle(label = names(genesets)[i]) + theme(plot.title = element_text(size = 12, face = "bold"))
  p2$labels["colour"] <- "signature mean"
  #ggsave(paste0(dir.name, "/", folders[5], "/", names(genesets)[i], "_clustree_mean.pdf"), plot = p2, scale = 1)
  
  #P-val plot
  p3 <- clustree(seurat, paste0(assay_type,"_snn_res."), node_colour = gs_name , node_colour_aggr = "mean") + ggtitle(label = names(genesets)[i]) + scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1.3) + theme(plot.title = element_text(size = 12, face = "bold"))
  clustree_gs_colname = paste0("mean_",gs_name)
  p3$data[[clustree_gs_colname]] <- -log(df_pval[, names(genesets)[i]])
  p3$labels[["colour"]] <- "-log(p_value)"
  #ggsave(paste0(dir.name, "/", folders[5], "/", names(genesets)[i], "_clustree_pval.pdf"), plot = p3, scale = 1) 
  
  #Plot combination
  p4 <- p2 + p3
  ggsave(paste0(dir.name, "/", folders[5], "/", names(genesets)[i], "_clustree_mean_pval.pdf"), plot = p4, scale = 1, width = 12, height=8, units = "in")
  
}

#9.9 Save seurat object
saveRDS(seurat, file = paste0(dir.name, "/",folders[5], "/seurat_complete.rds"))
