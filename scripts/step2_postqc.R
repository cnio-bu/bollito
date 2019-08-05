library("Seurat")
library("dplyr")
library("data.table")
library("reticulate")
library("ggplot2")

# A. Parameters: folder configuration 
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_celltypeid", "3_postprocessing", "4_degs", "5_gs")

# B. Parameters: analysis configuration 
# project_name = "Test"
min = as.numeric(snakemake@params[["min"]])
max = as.numeric(snakemake@params[["max"]])
mit = as.numeric(snakemake@params[["mit"]])
ribo = as.numeric(snakemake@params[["ribo"]])
filter.out = c(snakemake@params[["filter_out"]]) # Check this
filter.threshold = snakemake@params[["filter_threshold"]]# -1 would mean "<", the rest means ">"

# C. Analysis
# Read RDS file from previous step
seurat = readRDS(paste0(dir.name, "/", folders[1], "/seurat_pre-qc.rds"))
# 3. We should apply the filterings once the QC plots (GenePlot and Violin plots) have been checked.
#seurat <- subset(seurat, subset = nFeature_RNA > min & nFeature_RNA < max & percent.mt < mit & percent.ribo < ribo )
cells_seurat <- FetchData(object = seurat, vars = "nFeature_RNA")
seurat <- seurat[, which(x = cells_seurat > min & cells_seurat < max)] 
mit_seurat <- FetchData(object = seurat, vars = "percent.mt")
seurat <- seurat[, which(x = mit_seurat < mit)]
ribo_seurat <- FetchData(object = seurat, vars = "percent.ribo")
seurat <- seurat[, which(x = ribo_seurat < ribo)]

# 3.1. QC: violin plots - After
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.25) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[1], "/4_vlnplot_ngene_numi_pctmit_afterfilt.png"), scale = 1.5)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"), ncol = 3, pt.size = 0.25) + theme(legend.position="bottom") 
ggsave(paste0(dir.name, "/",folders[1], "/5_vlnplot_ngene_numi_pctribo_afterfilt.png"), scale = 1.5)

# 4. If there are negative markers availale: filter out cells based on gene expression. In this specific case, we are filtering out all cells expressing: Epcam, Pecam1, Krt19 and Ptprc. CHECK THIS

if(length(filter.out) > 0){
	if (filter.threshold == 0){
		for(i in 1:length(filter.out)){
			sub_seurat <- FetchData(object = seurat, vars = filter.out[i])
			seurat <- seurat[, which(x = sub_seurat > filter.threshold)]
		}			
	} else {
		for(i in 1:length(filter.out)){
			sub_seurat <- FetchData(object = seurat, vars = filter.out[i])
			seurat <- seurat[, which(x = sub_seurat == 0)]
		}	
	}
} else {
	next()
}

# Save RDS: we can use this object to generate all the rest of the data
saveRDS(seurat, file = paste0(dir.name, "/",folders[1], "/seurat_post-qc.rds"))

