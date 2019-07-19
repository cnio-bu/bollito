library("Seurat")
library("dplyr")
library("data.table")
library("reticulate")
library("ggplot2")
library("qusage")

# A. Parameters: folder configuration 
data_dir = paste0(snakemake@params[["input_dir"]],"/","Solo.out")
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_celltypeid", "3_postprocessing", "4_degs", "5_gs")

# B. Parameters: analysis configuration 
cell_cycle_file = snakemake@params[["cc_file"]]
geneset_collection = snakemake@params[["gs_collection"]]

# C. Analysis
seurat <- readRDS(paste0(dir.name, "/", folders[4], "/seurat_degs.rds"))
# 10 Cell cycle
dir.create(paste0(dir.name, "/", folders[5]))
# Read in a list of cell cycle markers, from Tirosh et al, 2015.
# We can segregate this list into markers of G2/M phase and markers of S phase.
cc.genes <- readLines(cell_cycle_file)
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

seurat <- CellCycleScoring(object = seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
FeaturePlot(object = seurat, features ="S.Score")
ggsave(paste0(dir.name, "/", folders[5], "/1_sscore_featureplot.pdf"))
FeaturePlot(object = seurat, features ="G2M.Score")
ggsave(paste0(dir.name, "/", folders[5], "/2_g2mscore_featureplot.pdf"))

# 11. Other genesets
genesets <- read.gmt(geneset_collection) #should be a tab file, each column = pathway.
seurat <- AddModuleScore(object = seurat, features= genesets, name = names(genesets))

for (i in 1:length(genesets)){
	module_name = colnames(seurat@meta.data)[grep(names(genesets)[i], colnames(seurat@meta.data))]
	FeaturePlot(object = seurat, features = module_name)
	ggsave(paste0(dir.name, "/", folders[5], "/", names(genesets)[i], "_featureplot.pdf"))
}

saveRDS(seurat, file = paste0(dir.name, "/",folders[5], "/seurat_complete.rds"))
