log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("Seurat"))

# A. Parameters: folder configuration.
dir.name = snakemake@params[["output_dir"]]
input_data = snakemake@input[["data"]]
random_seed = snakemake@params[["random_seed"]]
velocyto = snakemake@params[["velocyto"]]
outdir_config = snakemake@params[["outdir_config"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs", "6_traj_in", "7_func_analysis", "1_merged")

# C. Analysis.
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}

# 2.5. Merging all samples in the execution. 

# Loading seurat object function.
combine_object <- function(x, velocyto) {
  experiment = head(tail(strsplit(x, "/")[[1]], n=3), n=1) #Assay name is stored to later use in integrated object metadata.
  seurat_obj <- readRDS(x)
  seurat_obj <- AddMetaData(object = seurat_obj,
                            metadata = experiment,
                            col.name = 'assay_name')
  # If RNA velocity is going to be performed, we add the velocyto matrices in this step.
  if (velocyto){
    # Velocyto matrices path is infered. 
    velocyto_dir = paste0(outdir_config,"/star/", experiment,"/Solo.out/Velocyto/raw/")
    velo_names = c("spliced", "unspliced", "ambiguous")
    vel_matrices = list()
    # The matrices are read in 10x format.
    for (name in velo_names) {
      vel_matrices[[name]] <- Read10X(data.dir = paste0(velocyto_dir, name))
    }
    # The matrices are added as assays in the respective seurat object.
    for (name in velo_names) {
      vel_matrices[[name]] <- vel_matrices[[name]][, which(x = colnames(vel_matrices[[name]]) %in% colnames(seurat_obj))] 
      seurat_obj[[name]] <- CreateAssayObject(counts = vel_matrices[[name]])
    }
  }
  return(seurat_obj)
}

# Function to get sample names.
get_name_assays <- function(x) {
  experiment = head(tail(strsplit(x, "/")[[1]], n=3), n=1) #Assay name is stored to later use in integrated object metadata.
  return(experiment)
}

# Seurat objects are loaded and sample names are obtained.
seurat_list = lapply(input_data, function(x) combine_object(x, velocyto))
cells_id = sapply(input_data, function(x) get_name_assays(x))
names(seurat_list) <- cells_id

# Seurat object split for merging.
x_seurat <- seurat_list[[1]]
y_seurat <- seurat_list[2:length(seurat_list)]

# Merging step.
seurat <- merge(x_seurat, y = as.vector(y_seurat), add.cell.ids = cells_id, project = "merged")

# Save expression matrix.
write.table(as.matrix(seurat@assays$RNA@counts), file = paste0(dir.name, "/", folders[1], "/expression_matrix.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Save merged Seurat object.
saveRDS(seurat, file = paste0(dir.name, "/", folders[1], "/seurat_post-qc.rds"))



