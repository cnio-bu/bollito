library("Seurat")
library("dplyr")
library("reticulate")
py_install("umap-learn") # UMAP needs to be installed in order to run the RunUMAP function

# Read STARSolo output
data_dir = "Solo.out/"
expression_matrix <- Read10X(data.dir = data_dir)

# Seurat folder
sys_date = Sys.Date()
dir.create(paste0(data_dir, sys_date,"_Seurat_Analysis"))
dir.name = paste0(data_dir, sys_date ,"_Seurat_Analysis")
folders = c("1_Preprocessing", "2_CellTypeID", "3_Postprocessing", "4_DEGs", "5_Cell_cycle")

# A.1. Beginning with Seurat: http://satijalab.org/seurat/
## Creating a seurat object 
seurat = CreateSeuratObject(expression_matrix, project = "project-name", min.features = 200) # Here, a project name should be added by the user. 

# A.2. Preprocessing: Filter out low-quality cells
dir.create(paste0(dir.name, "/", folders[1]))

## 1. Mitochondrial genes - check levels of expression for mt genes 
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
## 2. Ribosomal genes - check levels of expression for rb genes 
seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^Rp[sl][[:digit:]]")

## 3. QC: violin plots
pdf(paste0(dir.name, "/",folders[1], "/VlnPlot_distr_nGene_nUMI_percentMit_percentRibo_beforeFiltering.pdf"), paper= "USr", width = 20)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.5)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"), ncol = 3, pt.size = 0.5)
dev.off()
## QC: GenePlot
pdf(paste0(dir.name, "/", folders[1], "/GenePlot_nUMI_vs_percentMit_nGene.pdf"), paper = "USr", width = 20)
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5)
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

## 4. We should apply these filterings once the QC plots (GenePlot and Violin plots) have been checked.
seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5 & percent.ribo < 40 )

## 6. QC: violin plots - After
pdf(paste0(dir.name, "/",folders[1], "/VlnPlot_distr_nGene_nUMI_percentMit_afterFiltering.pdf"), paper= "USr", width = 20)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.5)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"), ncol = 3, pt.size = 0.5)
dev.off()

## 7. If there are negative markers availale: filter out cells based on gene expression. In this specific case, we are filtering out all cells expressing: Epcam, Pecam1, Krt19 and Ptprc. 
seurat <- subset(seurat, subset = Epcam > 0 | Pecam1 > 0 | Krt19 > 0 | Ptprc > 0 , invert = TRUE) 

## 8. Normalize data
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)

## 9. Find Variable Genes
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2500)
## Identify the 10 most highly variable genes
pdf(paste0(dir.name, "/",folders[1], "/VariableFeaturesPlot.pdf"), paper= "USr", width = 20)
top10 <- head(VariableFeatures(seurat), 15)
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()


# A.3. Start of Identifying Cell Types
dir.create(paste0(dir.name, "/", folders[2]))

## 1. Scale the data
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)

## 2. Run PCA
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), npcs = 100) # This result could all be saved in a table. 
# print(seurat[["pca"]], dims = 1:5, nfeatures = 5)
## Visualizing PCA in Different Ways: elbow plot most variable genes 
pdf(paste0(dir.name, "/", folders[3], "/VizPCA_most_Vargenes.pdf"), paper = "USr", width = 14)
VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
DimPlot(seurat, reduction = "pca")
DimHeatmap(seurat, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

## 3. Determine the dimensionality of the dataset
seurat <- JackStraw(seurat, num.replicate = 100, dims = 100)
seurat <- ScoreJackStraw(seurat, dims = 1:100)
## Visualize
pdf(paste0(dir.name, "/", folders[3], "/Significant_PCs.png"))
JackStrawPlot(seurat, dims = 1:100)
ElbowPlot(seurat, ndims = 100)
dev.off()

# A.4. Post-processing
dir.create(paste0(dir.name, "/", folders[3]))

## 1. FindClusters
set.seed(8458)
seurat <- FindNeighbors(seurat, dims = 1:50)
seurat <- FindClusters(seurat)
head(Idents(seurat), 5)

## 2. Run UMAP
seurat <- RunUMAP(seurat, dims = 1:50)
## Visualize Umap
pdf(paste0(dir.name, "/", folders[3], "/UMAP_min.pdf"), paper = "USr", width = 12)
DimPlot(seurat, reduction = "umap",  pt.size = 0.75, label = TRUE, label.size = 5) + theme_minimal()
dev.off()
## Visualize no - Umap
seurat.no.umap <- seurat
seurat.no.umap[["umap"]] <- NULL
pdf(paste0(dir.name, "/", folders[3], "/noUMAP.pdf"), paper = "USr", width = 15)
DimPlot(seurat.no.umap, pt.size = 1, label = TRUE, label.size = 5) + RotatedAxis()
dev.off()

# A.5 Differentially expressed genes between clusters. 
dir.create(paste0(dir.name, "/", folders[4]))

## 1. Table on differentially expressed genes - using basic filterings
for (i in 1:length(unique(Idents(seurat)))){
  clusterX.markers <- FindMarkers(seurat, ident.1 = unique(Idents(seurat))[i], min.pct = 0.25) #min expressed
  print(x = head(x = clusterX.markers, n = 5))
  write.table(clusterX.markers, file=paste0(dir.name, "/", folders[4], "/cluster",unique(Idents(seurat))[i],".markers.txt"), sep="\t", col.names = NA)
}

## 2. DE includying all genes - needed for a GSEA analysis. 
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

## 3. Find TOP markers
seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top.markers = seurat.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
## Visualize
pdf(paste0(dir.name, "/", folders[4], "/Vln_levels_expr_rawCounts_normCounts.pdf"), width = 14, height = 10)
VlnPlot(seurat, features=top.markers, slot="counts", log=TRUE)
VlnPlot(seurat, features=top.markers)
dev.off()
# Visualize co-expression of two features simultaneously
pdf(paste0(dir.name, "/", folders[4], "/FeaturePlot_topMarkers.pdf"), width = 25 , height = 7)
FeaturePlot(seurat, features = top.markers, blend = TRUE, pt.size = 1) 
dev.off()


## 4. Heatmap Variable Features
pdf(paste0(dir.name, "/", folders[4], "/Heatmap_Variable_Features.pdf"), width = 16, height = 20)
DoHeatmap(seurat, features = VariableFeatures(seurat)[1:100], cells = 1:1000, size = 4, angle = 90)
dev.off()


saveRDS(seurat, file = paste0(dir.name, "/", "seurat_final.rds"))


