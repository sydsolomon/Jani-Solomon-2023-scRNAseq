library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)

dark2 <- c("#1B9E77", "#D95F02", "#E7298A", "#7570B3", "#E6AB02","#66A61E") # sample color scheme

# Load Data and create Seurat Object
# Load in UMI matrix
umis <- Read10X_h5("filtered_feature_bc_matrix_gex.h5")

# Load in HTO matrix
hto1 <- Read10X_h5("filtered_feature_bc_matrix_hto1.h5")

# Select cell barcodes detected by both RNA & HTO
joint.bcs <- intersect(colnames(umis), colnames(hto1))
length(joint.bcs)

#Subset RNA & HTO counts by joint cell barcodes
umis <- umis[, joint.bcs]
hto1 <- hto1[,joint.bcs]

# Confirm correct HTO names
rownames(hto1) <- c("WT_Unstim", "WT_PIM6", "WT_H37Rv", "TLR2_Unstim", "TLR2_PIM6", "TLR2_H37Rv")

# Setup Seurat Object
sobj <- CreateSeuratObject(umis)

# Apply sctransform normalization 
# store michondrial percentage in metadata
sobj <- PercentageFeatureSet(sobj, pattern = "^mt-", col.name = "percent.mt")

# run sctransform & regress on mitochondrial percentage
sobj <- SCTransform(sobj, vars.to.regress = "percent.mt", verbose = FALSE)

## PCA & UMAP
# Use more PCs with sctransform
sobj <- RunPCA(sobj, verbose = FALSE)
sobj <- RunUMAP(sobj, dims = 1:30, verbose = FALSE)

sobj <- FindNeighbors(sobj, dims = 1:30, verbose = FALSE)
sobj <- FindClusters(sobj, resolution = 0.6, verbose = FALSE) #Optimal cluster resolution after testing resolution = c(0.4, 0.6, 0.8, 1, 1.2)

Idents(sobj) <- "SCT_snn_res.0.6" #Best cluster resolution
DimPlot(sobj, reduction = "umap", label = TRUE, label.size = 6, pt.size = 1) #check cluster 6 for mitochondrial expression due to spread
cluster6.markers <- FindMarkers(sobj, ident.1 = 6, min.pct = 0.25)
head(cluster6.markers, 20) # mitochondrial driven - remove

# Add HTO information
sobj[["HTO"]] <- CreateAssayObject(counts = hto1)
sobj # confirm addition of HTO assay

# Normalize HTO counts
sobj <-  NormalizeData(sobj, assay = "HTO", normalization.method = "CLR")

# Sample demux/doublets
sobj <- HTODemux(sobj, positive.quantile = 0.97)

# Filter cells (mt reads, cluster 6, doublets)
filt_so <- subset(sobj,
                  (nFeature_RNA > 250) &
                    (percent.mt < 20))
Idents(filt_so) <- "SCT_snn_res.0.6" 
filt_so <- subset(filt_so, idents = "6", invert = TRUE)
Idents(filt_so) <- "HTO_classification.global"
filt_so <- subset(filt_so, idents = "Doublet", invert = TRUE)
Idents(filt_so) <- "HTO_classification"
levels(x = filt_so) <- c("WT-Unstim", "WT-PIM6", "WT-H37Rv", "TLR2-Unstim", "TLR2-PIM6", "TLR2-H37Rv")
my_levels <- c("WT-Unstim", "WT-PIM6", "WT-H37Rv", "TLR2-Unstim", "TLR2-PIM6", "TLR2-H37Rv")
filt_so@meta.data[["HTO_classification"]] <- factor (x = filt_so@meta.data[["HTO_classification"]], levels = my_levels)

# Confirm sample identity assignment
RidgePlot(filt_so, assay = "HTO", features = rownames(filt_so[["HTO"]])[1:6], ncol = 3, cols = dark2)

VlnPlot(
  filt_so,
  features = c("nCount_RNA", "nFeature_RNA"),
  ncol = 2,
  pt.size = FALSE,
  cols = dark2)

# UMAP Visualization - Fig 4A
DimPlot(filt_so, pt.size = 1, cols = dark2)
FeaturePlot(filt_so, features = 'hto_WT-Unstim', cols = c("lightgrey", dark2[1]), pt.size = 1)
FeaturePlot(filt_so, features = 'hto_WT-PIM6', cols = c("lightgrey", dark2[2]), pt.size = 1)
FeaturePlot(filt_so, features = 'hto_WT-H37Rv', cols = c("lightgrey", dark2[3]), pt.size = 1)
FeaturePlot(filt_so, features = 'hto_TLR2-Unstim', cols = c("lightgrey", dark2[4]), pt.size = 1)
FeaturePlot(filt_so, features = 'hto_TLR2-PIM6', cols = c("lightgrey", dark2[5]), pt.size = 1)
FeaturePlot(filt_so, features = 'hto_TLR2-H37Rv', cols = c("lightgrey", dark2[6]), pt.size = 1)

# Principal component analysis
DimPlot(filt_so, reduction = "pca", label = FALSE, cols = dark2, pt.size = 1) # Supp Fig S4A 

# PCA Loadings for Supplemental Table X
# write.csv(Loadings(filt_so, reduction = "pca")[,1:2], './file.csv')

# Eigenvalues for PC Variance
pca <- filt_so[["pca"]]
eigValues <-(pca@stdev)^2 ## Eigenvalues
varExplained <- eigValues / sum(eigValues) # Export csv with values and plot in Prism for Supp Fig S4B

# scPC Scoring top 50 positively contributing genes to scPC 1 and 2
# Extract loadings for PC genes
PC_loadings <- Loadings(filt_so[["pca"]] )
PC_loadings <- PC_loadings[,1:2] # Reduce size for only PC1 and PC2
PC1_pos <- subset(PC_loadings, PC_loadings[,'PC_1'] > 0)
PC1_pos_order <- PC1_pos[order(-PC1_pos[,'PC_1']),]
PC1_top50 <- PC1_pos_order[1:50,]

PC2_pos <- subset(PC_loadings, PC_loadings[,'PC_2'] > 0)
PC2_pos_order <- PC2_pos[order(-PC2_pos[,'PC_2']),]
PC2_top50 <- PC2_pos_order[1:50,]

filt_so <- AddModuleScore(filt_so, features = rownames(PC1_top50), nbin = 30, ctrl = 30, name = "scPC1_Score")
filt_so <- AddModuleScore(filt_so, features = rownames(PC2_top50), nbin = 30, ctrl = 30, name = "scPC2_Score")

VlnPlot(filt_so, features = "scPC1_Score1", pt.size = FALSE, col = dark2)+ theme(legend.position = 'none') # Fig 4B
VlnPlot(filt_so, features = "scPC2_Score1", pt.size = FALSE, col = dark2)+ theme(legend.position = 'none') # Fig 4C

# Expression of individual scPC1 & scPC2 genes
VlnPlot(filt_so, features = "Tnf", pt.size = FALSE, col = dark2)+ theme(legend.position = 'none') # Tnf Fig 4D
# scPC1
VlnPlot(filt_so, features = 'Il1b',col = dark2, pt.size = FALSE)+ theme(legend.position = 'none')
VlnPlot(filt_so, features = 'Saa3', col = dark2, pt.size = FALSE)+ theme(legend.position = 'none')
VlnPlot(filt_so, features = 'Cxcl2', col = dark2, pt.size = FALSE)+ theme(legend.position = 'none')
scPC1_exp <- DotPlot(filt_so, features = c("Il1b", "Saa3", "Cxcl2","Tnf"))+ scale_y_discrete(limits = rev(levels(filt_so$HTO_classification)))
scPC1_exp

#scPC2
VlnPlot(filt_so, features = 'Cxcl10',col = dark2, pt.size = FALSE)+ theme(legend.position = 'none')
VlnPlot(filt_so, features = 'Rsad2',col = dark2, pt.size = FALSE)+ theme(legend.position = 'none')
VlnPlot(filt_so, features = 'Isg15',col = dark2, pt.size = FALSE)+ theme(legend.position = 'none')
VlnPlot(filt_so, features = 'Ifit2', col = dark2, pt.size = FALSE) + theme(legend.position = 'none')
scPC2_exp <- DotPlot(filt_so, features = c('Cxcl10', 'Rsad2', 'Isg15', 'Ifit2'))+ scale_y_discrete(limits = rev(levels(filt_so$HTO_classification)))
scPC2_exp

# Total cells in each score
pc1_cells <-  which(FetchData(filt_so, vars = 'scPC1_Score1')>1)
pc2_cells <-  which(FetchData(filt_so, vars = 'scPC2_Score1')>1)

# Identify sample identity of cells
pc1_identity <- table(filt_so@meta.data$HTO_classification[pc1_cells])
pc2_identity <- table(filt_so@meta.data$HTO_classification[pc2_cells])

# Coexpression of scPC1 & scPC2 gene module
FeaturePlot(filt_so, features = c('scPC1_Score1', 'scPC2_Score1'), pt.size = 1, blend = TRUE) # Fig 4F

pc1_pc2_cells <- which(FetchData(filt_so, vars = 'scPC1_Score1')>1 & FetchData(filt_so, vars = 'scPC2_Score1')>1) #row of cells that meet criteria
pc1_pc2_identity <- table(filt_so$HTO_classification[pc1_pc2_cells])
pc1_pc2_cells <- rownames(filt_so@meta.data[pc1_pc2_cells,]) # convert row number to cell barcode
DimPlot(filt_so, cells.highlight = pc1_pc2_cells, cols.highlight = "#FF00FF", cols = 'gray', order = TRUE)+ theme(legend.position = 'none')

# Cluster Analysis - Figure 5
Idents(filt_so) <- "SCT_snn_res.0.6" # Figure 5A
DimPlot(filt_so, group.by = "SCT_snn_res.0.6", pt.size = 1, label = TRUE, label.size = 6)

# Determine percentage of cluster of each sample - Figure 5B
# Subset by HTO then determine how many cells in each cluster
Idents(filt_so) <- "hash.ID" 
WT_Unstim <- subset(filt_so, idents = "WT-Unstim")
WT_PIM6 <- subset(filt_so, idents = "WT-PIM6")
WT_H37Rv <- subset(filt_so, idents = "WT-H37Rv")
TLR2_Unstim <- subset(filt_so, idents = 'TLR2-Unstim')
TLR2_PIM6 <- subset(filt_so, idents = 'TLR2-PIM6')
TLR2_H37Rv <- subset(filt_so, idents = 'TLR2-H37Rv')

Idents(WT_Unstim) <- 'SCT_snn_res.0.6'
Idents(WT_PIM6) <- 'SCT_snn_res.0.6'
Idents(WT_H37Rv) <- 'SCT_snn_res.0.6'
Idents(TLR2_Unstim) <- 'SCT_snn_res.0.6'
Idents(TLR2_PIM6) <- 'SCT_snn_res.0.6'
Idents(TLR2_H37Rv) <- 'SCT_snn_res.0.6'
table(Idents(WT_Unstim))
table(Idents(WT_PIM6))
table(Idents(WT_H37Rv))
table(Idents(TLR2_Unstim))
table(Idents(TLR2_PIM6))
table(Idents(TLR2_H37Rv))

# Identify cluster markers - Fig 5C-D
Idents(filt_so) <- 'SCT_snn_res.0.6'
cluster0.markers <- FindMarkers(filt_so, ident.1 = "0", min.pct = 0.25)
cluster1.markers <- FindMarkers(filt_so, ident.1 = "1", min.pct = 0.25)
cluster2.markers <- FindMarkers(filt_so, ident.1 = "2", min.pct = 0.25)
cluster3.markers <- FindMarkers(filt_so, ident.1 = '3', min.pct = 0.25)
cluster4.markers <- FindMarkers(filt_so, ident.1 = '4', min.pct = 0.25)
cluster5.markers <- FindMarkers(filt_so, ident.1 = '5', min.pct = 0.25)
cluster7.markers <- FindMarkers(filt_so, ident.1 = '7', min.pct = 0.25) 
cluster8.markers <- FindMarkers(filt_so, ident.1 = '8', min.pct = 0.25)
# Export marker genes to csv

# Individual cluster 7 genes - Figure 4E-F
Idents(filt_so) <- "HTO_classification"
VlnPlot(filt_so, features = 'Fabp4',col = dark2, pt.size = FALSE)+ theme(legend.position = 'none') 
VlnPlot(filt_so, features = 'Fabp5',col = dark2, pt.size = FALSE)+ theme(legend.position = 'none')
VlnPlot(filt_so, features = 'Gpnmb',col = dark2, pt.size = FALSE)+ theme(legend.position = 'none')
VlnPlot(filt_so, features = 'Lgals4',col = dark2, pt.size = FALSE)+ theme(legend.position = 'none') 
VlnPlot(filt_so, features = 'Apoe',col = dark2, pt.size = FALSE)+ theme(legend.position = 'none')
VlnPlot(filt_so, features = 'Il7r',col = dark2, pt.size = FALSE)+ theme(legend.position = 'none')
scCluster7_exp <- DotPlot(filt_so, features = c("Fabp4", "Fabp5", "Gpnmb", "Lgals4", "Apoe","Il7r")) + scale_y_discrete(limits = rev(levels(filt_so$HTO_classification)))
scCluster7_exp

# Cluster 7 score
cluster7.fc <- cluster7.markers[order(-cluster7.markers[,'avg_log2FC']),]
cluster7.fc50 <- cluster7.fc[1:50,]
filt_so <- AddModuleScore(filt_so, features = rownames(cluster7.fc50), nbin = 30, ctrl = 30, name = "Cluster7_Score")

VlnPlot(filt_so, features = "Cluster7_Score1", pt.size = FALSE, col = dark2) +theme(legend.position = 'none')

# Coexpression of scPC1 & scPC2 with Cluster7 gene expressing cells - Supplemental Figure S5A-B
FeaturePlot(filt_so, features = c('scPC1_Score1', 'Cluster7_Score1'), pt.size = 1, blend = TRUE)
FeaturePlot(filt_so, features = c('scPC2_Score1', 'Cluster7_Score1'), pt.size = 1, blend = TRUE)

# Co-expression cell number
cluster7_cells <-  which(FetchData(filt_so, vars = 'Cluster7_Score1')>1) # number of cells in cluster 7
pc1_cluster7_cells <- which(FetchData(filt_so, vars = 'scPC1_Score1')>1 & FetchData(filt_so, vars = 'Cluster7_Score1')>1)
pc2_cluster7_cells <- which(FetchData(filt_so, vars = 'scPC2_Score1')>1 & FetchData(filt_so, vars = 'Cluster7_Score1')>1)
pc1_pc2_cluster7_cells <- which(FetchData(filt_so, vars = 'scPC1_Score1')>1 & FetchData(filt_so, vars = 'scPC2_Score1')>1 & 
                                  FetchData(filt_so, vars = 'Cluster7_Score1') > 1) # no cells express all 3 modules
# Identify sample identity of cells
cluster7_identity <- table(filt_so$HTO_classification[cluster7_cells])
pc1_cluster7_identity <- table(filt_so$HTO_classification[pc1_cluster7_cells])
pc2_cluster7_identity <- table(filt_so$HTO_classification[pc2_cluster7_cells])

# Convert row number to cell barcode
pc2_cluster7_cells <- rownames(filt_so@meta.data[pc2_cluster7_cells,])
pc1_cluster7_cells <- rownames(filt_so@meta.data[pc1_cluster7_cells,])
DimPlot(filt_so, cells.highlight = pc1_cluster7_cells, cols.highlight = "#FF00FF", cols = 'gray', order = TRUE) + theme(legend.position = 'none')
DimPlot(filt_so, cells.highlight = pc2_cluster7_cells, cols.highlight = "#FF00FF", cols = 'gray', order = TRUE) + theme(legend.position = 'none')
