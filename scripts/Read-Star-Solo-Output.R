library("Seurat")
library("dplyr")

# Read STARSolo output
data_dir = "Solo.out/"
expression_matrix <- Read10X(data.dir = data_dir)
# Seurat folder
dir.create(paste0(data_dir, Sys.Date(),"_Seurat_Analysis"))
dir.name = paste0(data_dir, Sys.Date(),"_Seurat_Analysis")
folders = c("1_Preprocessing", "2_CellTypeID", "3_Postprocessing", "4_DEGs", "5_Further_divisions", "6_Cell_cycle")

# A.1. Beginning with Seurat: http://satijalab.org/seurat/
## Creating a seurat object 
seurat = CreateSeuratObject(raw.data = expression_matrix)

# A.2. Preprocessing: Filter out low-quality cells
dir.create(paste0(dir.name, "/", folders[1]))
# Mitochondrial genes
mito.genes <- grep(pattern = "^mt-", x = rownames(x = seurat@data), value = TRUE)
percent.mito <- Matrix::colSums(seurat@raw.data[mito.genes, ])/Matrix::colSums(seurat@raw.data)
head(seurat@meta.data) # Before adding
# Metadata adding: MT genes. AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats.
seurat <- AddMetaData(object = seurat, metadata = percent.mito, col.name = "percent.mito")
head(seurat@meta.data) # After adding
# Ribosomal genes
ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = seurat@data), value = TRUE)
percent.ribo <- Matrix::colSums(seurat@raw.data[ribo.genes, ])/Matrix::colSums(seurat@raw.data)
seurat <- AddMetaData(object = seurat, metadata = percent.ribo, col.name = "percent.ribo")
pdf(paste0(dir.name, "/",folders[1], "/VlnPlot_distr_nGene_nUMI_percentMit_percentRibo_beforeFiltering.pdf"), paper= "USr", width = 14)
VlnPlot(object = seurat, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3,  point.size.use = 0.1)
VlnPlot(object = seurat, features.plot = c("nGene", "nUMI", "percent.ribo"), nCol = 3,  point.size.use = 0.1)
dev.off()

#GenePlot
pdf(paste0(dir.name, "/", folders[1], "/GenePlot_nUMI_vs_percentMit_nGene.pdf"), paper = "USr", width = 14)
par(mfrow = c(1, 2))
GenePlot(object = seurat, gene1 = "nUMI", gene2 = "percent.mito", cex.use = 0.5)
GenePlot(object = seurat, gene1 = "nUMI", gene2 = "nGene", cex.use = 0.5)
dev.off()

# Apply filterings
seurat <- FilterCells(object = seurat, subset.names = c("nGene","percent.mito", "percent.ribo"), low.thresholds = c(200, -Inf, -Inf), high.thresholds = c(4000, 0.15, 0.4))
pdf(paste0(dir.name, "/",folders[1], "/VlnPlot_distr_nGene_nUMI_percentMit_afterFiltering.pdf"), paper= "USr", width = 14)
VlnPlot(object = seurat, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3,  point.size.use = 0.1) 
VlnPlot(object = seurat, features.plot = c("nGene", "nUMI", "percent.ribo"), nCol = 3,  point.size.use = 0.1)
dev.off()

# Normalize data
seurat  <- NormalizeData(object = seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Find Variable Genes
seurat <- FindVariableGenes(object = seurat, 
                            mean.function = ExpMean, 
                            dispersion.function = LogVMR, 
                            do.plot = FALSE)


# A.3. Start of Identifying Cell Types
dir.create(paste0(dir.name, "/", folders[2]))

# Scaling
seurat <- ScaleData(object = seurat) #you can give the scaling variables to regress out: delete those variables from the model
View(seurat@meta.data) # All these variables can be regressed out. If you have a variable that has indicated the batch your data belongs to, you can also substract it, but we'll look into that during the batch correction module.

# Running PCA
# This will run pca on the 
seurat <- RunPCA(object = seurat, 
                 pc.genes = seurat@var.genes,  #PCA only on the variable genes, you could use a set of genes or pathway.
                 do.print = TRUE, 
                 pcs.print = 1:5,  # Print these variable genes out
                 genes.print = 5)

# Running ICA
seurat <- RunICA(seurat, ics.compute=5) # check pbmc@dr

# Visualizing PCA in Different Ways: elbow plot most variable genes 
pdf(paste0(dir.name, "/", folders[2], "/VizPCA_most_Vargenes.pdf"), paper = "USr", width = 14)
VizPCA(object = seurat, pcs.use = 1:2) 
dev.off()

# Visualizing ICA in Different Ways
pdf(paste0(dir.name, "/", folders[2], "/VizICA_most_Vargenes.pdf"), paper = "USr", width = 14)
VizICA(object = seurat, ics.use=1:3)
dev.off()
#pdf(paste0(dir.name, "/", folders[2], "/PCAPlot_byBatch.pdf"), paper = "USr", width = 14)
#PCAPlot(object = seurat, dim.1 = 1, dim.2 = 2)
#dev.off()

# A.4. Post-processing
dir.create(paste0(dir.name, "/", folders[3]))

# Genes by PCs
seurat <- RunPCA(seurat, pc.genes = seurat@var.genes, pcs.compute = 30, do.print = TRUE, pcs.print = 5, genes.print = 5)
seurat <- ProjectPCA(object = seurat, do.print = FALSE)

pdf(paste0(dir.name, "/", folders[3], "/FeaturePlot_pca.pdf"), width = 24, height = 14)
FeaturePlot(seurat, 'nGene', cols.use=c("grey", "blue"), pt.size = 3, no.legend = FALSE, reduction.use = "pca")
dev.off()

# DO NOT RUN! Heatmap of the PC and the differentially expressed genes within that PC
#PCHeatmap(object = seurat, 
          #pc.use = 1, 
          #do.balanced = TRUE, 
          #label.columns = FALSE)

# PCElbowPlot
pdf(paste0(dir.name, "/", folders[3], "/PCElbowPlot.pdf"), width = 24, height = 14)
PCElbowPlot(object = seurat,  num.pc = 30) # How to choose the number of PCs? Use this plot!! Maybe cut around 10, because is flattening from that on...? You can go even to 10-20.
dev.off()

# FindClusters
set.seed(2019)
seurat <- FindClusters(object = seurat, 
                          reduction.type = "pca", 
                          dims.use = 1:20, 
                          resolution = 0.6, 
                          print.output = 0) # TSNEPlot: 9 clusters found.
# for res = 1; fcnio.CCA <-RunTSNE(fcnio.CCA, reduction.use = "pca", dims.use = 1:10, perplexity=10) #Takes the first pcs, for example 1 to 5.

#tSNE
seurat <-RunTSNE(seurat, reduction.use = "pca", dims.use = 1:20, perplexity=10) #Takes the first pcs, for example 1 to 5.

pdf(paste0(dir.name, "/", folders[3], "/tSNE_FeaturePlot_Dimplot.pdf"), width = 20, height = 14)
TSNEPlot(object = seurat, pt.size =3, do.label=T, label.size = 10)
FeaturePlot(seurat, 'nGene', cols.use=c("grey", "blue"), pt.size = 3, no.legend = FALSE, reduction.use = "tsne")
#DimPlot(seurat, reduction.use = "tsne", dim.1 = 1, dim.2 = 2, group.by = "res.0.6", do.return = T, pt.size = 3)
dev.off()

# A.5 Differentially expressed genes between clusters. 
dir.create(paste0(dir.name, "/", folders[4]))

for (i in 1:length(unique(seurat@meta.data$res.0.6))){
  clusterX.markers <- FindMarkers(object = seurat, ident.1 = unique(seurat@meta.data$res.0.6)[i], min.pct = 0.25)
  print(x = head(x = clusterX.markers, n = 5))
  write.table(clusterX.markers, file=paste0(dir.name, "/", folders[4], "/cluster",unique(seurat@meta.data$res.0.6)[i],".markers.txt"), sep="\t", col.names = NA)
}

# Find markers for every cluster compared to all remainin cells. 
seurat.markers <- FindAllMarkers(object = seurat, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
groupedby.clusters.markers = seurat.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(groupedby.clusters.markers, file=paste0(dir.name, "/", folders[4], "/groupedby_clusters_markers.txt"), sep="\t", col.names = NA)
#top10 = seurat.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
#write.table(top10, file=paste0(dir.name, "/", folders[4], "/groupedby_clusters_markers_10.txt"), sep="\t", col.names = NA)

# Feature plot
pdf(paste0(dir.name, "/", folders[4], "/FeaturePlots_markers.pdf"), width = 20, height = 14)
TSNEPlot(object = seurat, pt.size = 3, do.label=T, label.size = 10)
FeaturePlot(object = seurat, features.plot = "Saa3", cols.use = c("grey", "blue"), pt.size = 3)
FeaturePlot(object = seurat, features.plot = "Lum", cols.use = c("grey", "blue"), pt.size = 3)
FeaturePlot(object = seurat, features.plot = "Hp", cols.use = c("grey", "blue"), pt.size = 3)
FeaturePlot(object = seurat, features.plot = "Pdgfra", cols.use = c("grey", "blue"), pt.size = 3)
FeaturePlot(object = seurat, features.plot = "Has1", cols.use = c("grey", "blue"), pt.size = 3)
FeaturePlot(object = seurat, features.plot = "Has2", cols.use = c("grey", "blue"), pt.size = 3)
FeaturePlot(object = seurat, features.plot = "Lrg1", cols.use = c("grey", "blue"), pt.size = 3)
FeaturePlot(object = seurat, features.plot = "Il6", cols.use = c("grey", "blue"), pt.size = 3)
#FeaturePlot(object = seurat, features.plot = "Sma", cols.use = c("grey", "blue"), pt.size = 3) NOT FOUND
FeaturePlot(object = seurat, features.plot = "Vim", cols.use = c("grey", "blue"), pt.size = 3)
FeaturePlot(object = seurat, features.plot = "Msln", cols.use = c("grey", "blue"), pt.size = 3)
dev.off()

# From https://www.nature.com/articles/s41467-018-07582-3/figures/1 NOT SAVED
#pdf(paste0(dir.name, "/", folders[4], "/FeaturePlots_markers_Nature_Bartoschek.pdf"), width = 20, height = 14)
TSNEPlot(object = seurat, pt.size = 3, do.label=T, label.size = 10)
FeaturePlot(object = seurat, features.plot = "Fap", cols.use = c("grey", "blue"), pt.size = 3)
FeaturePlot(object = seurat, features.plot = "S100a4", cols.use = c("grey", "blue"), pt.size = 3)
FeaturePlot(object = seurat, features.plot = "Sparc", cols.use = c("grey", "blue"), pt.size = 3)
FeaturePlot(object = seurat, features.plot = "Acta2", cols.use = c("grey", "blue"), pt.size = 3)
FeaturePlot(object = seurat, features.plot = "Pdgfra", cols.use = c("grey", "blue"), pt.size = 3)
FeaturePlot(object = seurat, features.plot = "Pdgfrb", cols.use = c("grey", "blue"), pt.size = 3)
FeaturePlot(object = seurat, features.plot = "Tnc", cols.use = c("grey", "blue"), pt.size = 3)
#dev.off()

require(lattice)
for (i in 1:length(unique(seurat@meta.data$res.0.6))){
  pdf(paste0(dir.name, "/", folders[4], "/VlnPlots_cluster", unique(seurat@meta.data$res.0.6)[i],"markers.pdf"), width = 24, height = 14)
  print(VlnPlot(object = seurat, features.plot = groupedby.clusters.markers$gene[groupedby.clusters.markers$cluster == unique(seurat@meta.data$res.0.6)[i]]))
  dev.off()
}

# Joint FeaturesPlot
pdf(paste0(dir.name, "/", folders[4], "/FeaturePlot_All_Cluster_markers.pdf"), width = 20, height = 14)
FeaturePlot(object = seurat, 
            features.plot = groupedby.clusters.markers$gene,
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
dev.off()

pdf(paste0(dir.name, "/", folders[4], "/FeaturePlot_All_Cluster_Selectedmarkers.pdf"), width = 20, height = 14)
FeaturePlot(object = seurat, 
            features.plot = c("Saa3", "Lum", "Hp", "Pdgfra", "Has1","Has2", "Lrg1", "Il6", "Vim", "Msln"),
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
dev.off()

# HeatMap top10
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
pdf(paste0(dir.name, "/", folders[4], "/Heatmap_markers.pdf"), width = 24, height = 14)
DoHeatmap(object = seurat, genes.use = groupedby.clusters.markers$gene, 
          remove.key = TRUE, cex.col=1, cex.row =14, slim.col.label =TRUE, title ="Cluster marker genes")
dev.off()

# Further subdivisions within cell types. If you perturb some of our parameter choices above (for example, setting resolution=0.8 or changing the number of PCs), you might see one cluster subdivide into two groups. 
#You can explore this subdivision to find markers separating the two subsets. However, before reclustering (which will overwrite object@ident), we can stash our renamed identities to be easily recovered later.
dir.create(paste0(dir.name, "/", folders[5]))
# First lets stash our identities for later
seurat <- StashIdent(object = seurat, save.name = "ClusterNames_0.6")

seurat <- FindClusters(object = seurat, 
                       reduction.type = "pca", 
                       dims.use = 1:20, 
                       resolution = 1) 

pdf(paste0(dir.name, "/", folders[5], "/tsne_markers_res1.pdf"), width = 20, height = 14)
plot <- TSNEPlot(object = seurat, 
                 do.return = TRUE, 
                 no.legend = TRUE, 
                 do.label = TRUE, 
                 pt.size = 2.5,
                 label.size = 10)

plot
dev.off()

# DEGS
for (i in 1:length(unique(seurat@meta.data$res.1))){
  clusterX.markers <- FindMarkers(object = seurat, ident.1 = unique(seurat@meta.data$res.1)[i], min.pct = 0.25)
  print(x = head(x = clusterX.markers, n = 5))
  write.table(clusterX.markers, file=paste0(dir.name, "/", folders[5], "/cluster",unique(seurat@meta.data$res.1)[i],".markers.txt"), sep="\t", col.names = NA)
}

# Find markers
seurat.markers <- FindAllMarkers(object = seurat, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
groupedby.clusters.markers = seurat.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(groupedby.clusters.markers, file=paste0(dir.name, "/", folders[5], "/groupedby_clusters_markers.txt"), sep="\t", col.names = NA)
# Markers Vln 
for (i in 1:length(unique(seurat@meta.data$res.1))){
  pdf(paste0(dir.name, "/", folders[5], "/VlnPlots_cluster", unique(seurat@meta.data$res.1)[i],"markers.pdf"), width = 24, height = 14)
  print(VlnPlot(object = seurat, features.plot = groupedby.clusters.markers$gene[groupedby.clusters.markers$cluster == unique(seurat@meta.data$res.1)[i]]))
  dev.off()
}

## Visualize the expression of the first 5 marker genes on tSNE across the different batches using FeatureHeatmap.
#seurat <- SetAllIdent(seurat, id = "res.0.6")
#pdf(paste0(dir.name, "/", folders[3], "/Batch_cor_FeatureHeatmap_markers.pdf"), width = 24, height = 14)
#FeatureHeatmap(seurat, features.plot = rownames(seurat.markers)[1:5], pt.size = 2.5, key.position = "top", max.exp = 3)
#dev.off()


# Joint FeaturesPlot
pdf(paste0(dir.name, "/", folders[5], "/FeaturePlot_All_Cluster_markers.pdf"), width = 20, height = 14)
FeaturePlot(object = seurat, 
            features.plot = groupedby.clusters.markers$gene,
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
dev.off()

# HeatMap top10
pdf(paste0(dir.name, "/", folders[5], "/Heatmap_markers.pdf"), width = 24, height = 14)
DoHeatmap(object = seurat, genes.use = groupedby.clusters.markers$gene, 
          remove.key = TRUE, cex.col=1, cex.row =14, slim.col.label =TRUE, title ="Cluster marker genes")
dev.off()

# A.6. Cell cycle
dir.create(paste0(dir.name, "/", folders[6]))
seurat <- StashIdent(object = seurat, save.name = "ClusterNames_1")

# Read in a list of cell cycle markers, from Tirosh et al, 2015.
# We can segregate this list into markers of G2/M phase and markers of S phase.
cc.genes <- readLines(paste0(data_dir, "regev_lab_cell_cycle_genes.txt"))
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

seurat <- CellCycleScoring(object = seurat, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = T)
pdf(paste0(dir.name, "/", folders[6], "/FeaturePlot_cellcycleSignatures.pdf"), width = 20, height = 14)
FeaturePlot(object = seurat, features.plot ="S.Score", pt.size = 2.5)
FeaturePlot(object = seurat, features.plot ="G2M.Score", pt.size = 2.5)
dev.off()


## Selected pathways

#REACTOME_INTERFERON_SIGNALING = c("AAAS", "ADAR", "ARIH1", "B2M", "CAMK2A", "CAMK2B", "CAMK2D", "CD44", "CIITA", "DDX58", "EGR1", "EIF2AK2", "EIF4A1", "EIF4A2", "EIF4A3", "EIF4E", "EIF4E2", "EIF4E3", "EIF4G1", "EIF4G2", "EIF4G3", "FCGR1A", "FCGR1B", "FLNB", "GBP1", "GBP2", "GBP4", "GBP5", "GBP6", "GBP7", "HERC5", "HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DRB1", "HLA-DRB3", "HLA-DRB5", "HLA-F", "HLA-G", "HLA-K", "ICAM1", "IFI27", "IFI35", "IFI6", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3", "IFNA1", "IFNA10", "IFNA14", "IFNA16", "IFNA17", "IFNA2", "IFNA21", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNAR1", "IFNAR2", "IFNB1", "IFNG", "IFNGR1", "IFNGR2", "IP6K2", "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "JAK1", "JAK2", "KPNA1", "KPNA2", "KPNA3", "KPNA4", "KPNA5", "KPNB1", "LOC144383", "LOC441019", "LOC646981", "MAPK3", "MT2A", "MX1", "MX2", "NCAM1", "NEDD4", "NUP107", "NUP133", "NUP153", "NUP155", "NUP188", "NUP205", "NUP210", "NUP214", "NUP35", "NUP37", "NUP43", "NUP50", "NUP54", "NUP62", "NUP85", "NUP88", "NUP93", "NUPL1", "NUPL2", "OAS1", "OAS2", "OAS3", "OASL", "PIAS1", "PIN1", "PLCG1", "PML", "POM121", "PPM1B", "PRKCD", "PSMB8", "PTAFR", "PTPN1", "PTPN2", "PTPN6", "RAE1", "RANBP2", "RNASEL", "RPS27A", "RPS27AP11", "SEH1L", "SOCS1", "SOCS3", "SP100", "STAT1", "STAT2", "SUMO1", "TPR", "TRIM25", "TYK2", "UBA52", "UBA7", "UBE2E1", "UBE2L6", "UBE2N", "USP18", "VCAM1", "XAF1")
#REACTOME_INTERFERON_GAMMA_SIGNALING = c("B2M", "CAMK2A", "CAMK2B", "CAMK2D", "CD44", "CIITA", "FCGR1A", "FCGR1B", "GBP1", "GBP2", "GBP4", "GBP5", "GBP6", "GBP7", "HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DRB1", "HLA-DRB3", "HLA-DRB5", "HLA-F", "HLA-G", "HLA-K", "ICAM1", "IFNG", "IFNGR1", "IFNGR2", "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", "IRF7", "IRF8", "IRF9", "JAK1", "JAK2", "LOC441019", "LOC646981", "MT2A", "NCAM1", "OAS1", "OAS2", "OAS3", "OASL", "PIAS1", "PML", "PRKCD", "PTAFR", "PTPN1", "PTPN2", "PTPN6", "SOCS1", "SOCS3", "SP100", "STAT1", "SUMO1", "VCAM1")
#REACTOME_INTERFERON_ALPHA_BETA_SIGNALING = c("ADAR", "EGR1", "GBP2", "HLA-A", "HLA-B", "HLA-C", "HLA-F", "HLA-G", "HLA-K", "IFI27", "IFI35", "IFI6", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3", "IFNA1", "IFNA10", "IFNA14", "IFNA16", "IFNA17", "IFNA2", "IFNA21", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNAR1", "IFNAR2", "IFNB1", "IP6K2", "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "JAK1", "LOC144383", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "PSMB8", "PTPN1", "PTPN6", "RNASEL", "SOCS1", "SOCS3", "STAT1", "STAT2", "TYK2", "USP18", "XAF1")
#KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY = c("ATG12", "ATG5", "AZI2", "CASP10", "CASP8", "CHUK", "CXCL10", "CYLD", "DAK", "DDX3X", "DDX3Y", "DDX58", "DHX58", "FADD", "IFIH1", "IFNA1", "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA2", "IFNA21", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNB1", "IFNE", "IFNK", "IFNW1", "IKBKB", "IKBKE", "IKBKG", "IL12A", "IL12B", "IL8", "IRF3", "IRF7", "ISG15", "MAP3K1", "MAP3K7", "MAPK10", "MAPK11", "MAPK12", "MAPK13", "MAPK14", "MAPK8", "MAPK9", "MAVS", "NFKB1", "NFKBIA", "NFKBIB", "NLRX1", "OTUD5", "PIN1", "RELA", "RIPK1", "RNF125", "SIKE1", "TANK", "TBK1", "TBKBP1", "TMEM173", "TNF", "TRADD", "TRAF2", "TRAF3", "TRAF6", "TRIM25")
#KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY = c("ARAF", "BID", "BRAF", "CASP3", "CD244", "CD247", "CD48", "CHP", "CHP2", "CSF2", "FAS", "FASLG", "FCER1G", "FCGR3A", "FCGR3B", "FYN", "GRB2", "GZMB", "HCST", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-G", "HRAS", "ICAM1", "ICAM2", "IFNA1", "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA2", "IFNA21", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNAR1", "IFNAR2", "IFNB1", "IFNG", "IFNGR1", "IFNGR2", "ITGAL", "ITGB2", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DS1", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KLRC1", "KLRC2", "KLRC3", "KLRD1", "KLRK1", "KRAS", "LAT", "LCK", "LCP2", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "MICA", "MICB", "NCR1", "NCR2", "NCR3", "NFAT5", "NFATC1", "NFATC2", "NFATC3", "NFATC4", "NRAS", "PAK1", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PIK3R5", "PLCG1", "PLCG2", "PPP3CA", "PPP3CB", "PPP3CC", "PPP3R1", "PPP3R2", "PRF1", "PRKCA", "PRKCB", "PRKCG", "PTK2B", "PTPN11", "PTPN6", "RAC1", "RAC2", "RAC3", "RAET1E", "RAET1G", "RAET1L", "RAF1", "SH2D1A", "SH2D1B", "SH3BP2", "SHC1", "SHC2", "SHC3", "SHC4", "SOS1", "SOS2", "SYK", "TNF", "TNFRSF10A", "TNFRSF10B", "TNFRSF10C", "TNFRSF10D", "TNFSF10", "TYROBP", "ULBP1", "ULBP2", "ULBP3", "VAV1", "VAV2", "VAV3", "ZAP70")
#KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION = c("B2M", "CALR", "CANX", "CD4", "CD74", "CD8A", "CD8B", "CIITA", "CREB1", "CTSB", "CTSL1", "CTSS", "HLA-A", "HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-E", "HLA-F", "HLA-G", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA1L", "HSPA2", "HSPA4", "HSPA5", "HSPA6", "HSPA8", "IFI30", "IFNA1", "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA2", "IFNA21", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DS1", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KLRC1", "KLRC2", "KLRC3", "KLRC4", "KLRD1", "LGMN", "LOC652614", "LTA", "NFYA", "NFYB", "NFYC", "PDIA3", "PSME1", "PSME2", "PSME3", "RFX5", "RFXANK", "RFXAP", "TAP1", "TAP2", "TAPBP")
#fcnio.CCA <- AddModuleScore(fcnio.CCA, genes.list = list(REACTOME_INTERFERON_SIGNALING, REACTOME_INTERFERON_GAMMA_SIGNALING, REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY, KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY, KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION), ctrl.size = 20, enrich.name = c("REACTOME_INTERFERON_SIGNALING", "REACTOME_INTERFERON_GAMMA_SIGNALING", "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING", "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY", "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY", "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION"))

#pdf(paste0(dir.name, "/", folders[4], "/FeaturePlot_signatures.pdf"), width = 14, height = 8)
#FeaturePlot(object = fcnio.CCA, features.plot = "REACTOME_INTERFERON_SIGNALING1", pt.size = 2.5)
#FeaturePlot(object = fcnio.CCA, features.plot = "REACTOME_INTERFERON_GAMMA_SIGNALING2", pt.size = 2.5)
#FeaturePlot(object = fcnio.CCA, features.plot = "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING3", pt.size = 2.5)
#FeaturePlot(object = fcnio.CCA, features.plot = "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY4", pt.size = 2.5)
#FeaturePlot(object = fcnio.CCA, features.plot = "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY5", pt.size = 2.5)
#FeaturePlot(object = fcnio.CCA, features.plot = "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION6", pt.size = 2.5)
#dev.off()

## Save current progress.
save(seurat, file = paste0(dir.name, "/",Sys.Date(), "_seurat.Rda"))
# To load the data, run the following command.
load(paste0(dir.name, "/",Sys.Date(), "_seurat.Rda"))
