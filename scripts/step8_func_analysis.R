log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("Seurat"))
suppressMessages(library("ggplot2"))
suppressMessages(library("viridis"))
suppressMessages(library("devtools"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("VISION"))
suppressMessages(library("scales"))

dir.name = snakemake@params[["output_dir"]]
input_data = snakemake@input[["seurat_obj"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis")

# B. Parameters: analysis configuration
mol_signatures = snakemake@params[["mol_signatures"]]
meta_columns = snakemake@params[["meta_columns"]]
n_cores =  snakemake@params[["n_cores"]]
selected_res = snakemake@params[["selected_res"]]
random_seed = snakemake@params[["random_seed"]]
regress_out = snakemake@params[["regress_out"]]
vars_to_regress = snakemake@params[["vars_to_regress"]]
regress_cell_cycle = snakemake@params[["regress_cell_cycle"]]

# C. Analysis
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}
# 11. Vision functional analysis
# 11.1 Seurat object with clustering is loaded.
seurat <- readRDS(input_data)
assay_type <- seurat@active.assay
cluster_res <- paste0(assay_type, "_snn_res.", selected_res)
if (!(cluster_res %in% colnames(seurat@meta.data))){
  stop("Specified resolution is not available.")
}
meta_columns <-  c(meta_columns, cluster_res)

# 11.2 Vision object is created. 
# If seurat object is not a integration
message(colnames(seurat@meta.data))
if (seurat@active.assay != "integrated"){
# To avoid negative values with SCT normalization we use the standard, keeping the information from the previous analysis.
  if (seurat@active.assay == "SCT") {
    seurat@active.assay <- "RNA"
    seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
    all.genes <- rownames(seurat)
    if(regress_out == TRUE){
        if (regress_cell_cycle) {
            message ("1")
            seurat <- ScaleData(seurat, features = rownames(seurat), vars.to.regress = c(vars_to_regress, "S.Score", "G2M.Score"))
        } else {
            message("2")
            seurat <- ScaleData(seurat, features = rownames(seurat), vars.to.regress = vars_to_regress)
        }
    } else {
        seurat <- ScaleData(seurat, features = rownames(seurat))
        message("3")
    }
  }
  suppressMessages(vis <- Vision(seurat,
		  signatures = mol_signatures,
                  meta = seurat@meta.data[,meta_columns],
                  projection_methods = NULL))
} else {
  suppressMessages(vis <- Vision(rescale(seurat@assays$integrated@scale.data, c(0,10)),
                    signatures = mol_signatures,
                    meta = seurat@meta.data[,c(meta_columns, "assay_name")],
                    projection_methods = "UMAP"))
}

# 11.3 Vision analysis step
options(mc.cores = n_cores)
vis <- suppressMessages(analyze(vis))

# 11.4 Outputs generation
# 11.4.1 Projection values stored.
if (seurat@active.assay != "integrated"){
    projections <- vis@Projections$Seurat_umap
} else {
    projections <- vis@Projections[[1]]
}

# 11.4.2 Cluster plot in UMAP projection
getPalette <- colorRampPalette(brewer.pal(9,'Set1'))

p1 <- ggplot() + aes(x=projections[, 1], y=projections[, 2], color=vis@metaData[,cluster_res]) + geom_point(alpha=0.5) + xlab("UMAP_1") + ylab("UMAP_2") + ggtitle("Clusters representation UMAP") + labs(color='clusters') + scale_color_manual(values = getPalette(nlevels(vis@metaData[,cluster_res])))
suppressMessages(ggsave(paste0(dir.name, "/", folders[7], "/CLUSTERS_REPRESENTATION_UMAP_set1.pdf"), plot = p1))

# 11.4.3 Integration plot in UMAP projection
if (seurat@active.assay == "integrated"){
  p2 <- ggplot() + aes(x=projections[, 1], y=projections[, 2], color=vis@metaData[,"assay_name"]) + geom_point(alpha=0.5) + xlab("UMAP_1") + ylab("UMAP_2") + ggtitle("Integration representation UMAP") + labs(color='Assay') + scale_color_manual(values = getPalette(nlevels(vis@metaData[,"assay_name"])))
  suppressMessages(ggsave(paste0(dir.name, "/", folders[7], "/INTEGRATION_REPRESENTATION_UMAP_set1.pdf"), plot = p2))
}

# 11.4.4 Clusters vs Molecular Signatures statistics values 
for (i in 1:nlevels(vis@metaData[,cluster_res])){
  write.table(vis@ClusterComparisons$Signatures[[cluster_res]][[i]], file = paste0(dir.name, "/", folders[7], "/vision_table_res_", selected_res, "_cluster_", (i-1), ".tsv"), col.names = NA, quote = FALSE, sep = "\t")
}

# 11.4.5 Molecular Signatures statistics values
stats_DF <- data.frame("Consistency" = vis@LocalAutocorrelation$Signatures$C,
                       "p-values" =  vis@LocalAutocorrelation$Signatures$pValue, "FDR" = vis@LocalAutocorrelation$Signatures$FDR)
colnames(stats_DF) <- c("Consistency", "p-values", "FDR")
rownames(stats_DF) <- rownames(vis@LocalAutocorrelation$Signatures) 
write.table(stats_DF, file = paste0(dir.name, "/", folders[7], "/vision_param_values_per_geneset.txt"), col.names = NA, quote = FALSE, sep = "\t" )

# 11.4.6 Molecular signatures scores in UMAP projection.
for (genesig in colnames(vis@SigScores)){
  if (stats_DF[genesig,"FDR"] < 0.05) {
    p3 <- ggplot() + aes(x=projections[, 1], y=projections[, 2],  color=vis@SigScores[, genesig]) + geom_point() + ggtitle(genesig) +   scale_color_viridis(option = "D") + xlab("UMAP_1") + ylab("UMAP_2") + ggtitle(paste0(genesig, " - UMAP")) +  labs(color='vision score')
    suppressMessages(ggsave(paste0(dir.name, "/", folders[7], "/", genesig, "_UMAP.pdf"), plot = p3))
  }
}

# 11.5 Save Vision Object
saveRDS(vis, file = paste0(dir.name, "/", folders[7], "/vision_object.rds"))
