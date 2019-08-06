suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))
suppressMessages(library("qusage"))

# A. Parameters: folder configuration 
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_degs", "5_gs")

# B. Parameters: analysis configuration 
geneset_collection = snakemake@params[["gs_collection"]]

# C. Analysis
seurat <- readRDS(paste0(dir.name, "/", folders[4], "/seurat_degs.rds"))

dir.create(paste0(dir.name, "/", folders[5]))

# 10. GS scoring
genesets <- read.gmt(geneset_collection) #should be a tab file, each column = pathway.
seurat <- AddModuleScore(object = seurat, features= genesets, name = names(genesets))

for (i in 1:length(genesets)){
	module_name = colnames(seurat@meta.data)[grep(names(genesets)[i], colnames(seurat@meta.data))]
	FeaturePlot(object = seurat, features = module_name) #+ theme(legend.position="bottom") 
	ggsave(paste0(dir.name, "/", folders[5], "/", names(genesets)[i], "_featureplot.pdf"), scale = 1.5)
}

saveRDS(seurat, file = paste0(dir.name, "/",folders[5], "/seurat_complete.rds"))
