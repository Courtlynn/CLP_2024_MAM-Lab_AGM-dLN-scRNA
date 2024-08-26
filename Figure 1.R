# Creating Seurat object
# QC data
# Merge objects 
# Join layers 

rm(list=ls())
set.seed(1234)
# Load the necessary libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(Matrix)
library(dplyr)
library(readr)
library(glmGamPoi)

# Set working directory to the location of Seurat objects
setwd("~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Objects")

# Load the Seurat objects
Flag_Seurat <- readRDS("Flag_Seurat.RDS")
M229_Seurat <- readRDS("M229_Seurat.RDS")
PBS_Seurat <- readRDS("PBS_Seurat.RDS")
R848_Seurat <- readRDS("R848_Seurat.RDS")

# Function for performing QC and filtering
perform_qc <- function(seurat_obj) {
  # Calculate the percentage of mitochondrial genes
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Visualize QC metrics
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  # Filter cells based on QC metrics
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 4)
  
  return(seurat_obj)
}

# Perform QC and filtering on each Seurat object and save new object after QC
Flag_Seurat <- perform_qc(Flag_Seurat)
saveRDS(Flag_Seurat, file = ("Flag_Seurat_QC.RDS"))

M229_Seurat <- perform_qc(M229_Seurat)
saveRDS(M229_Seurat, file = ("M229_Seurat_QC.RDS"))

PBS_Seurat <- perform_qc(PBS_Seurat)
saveRDS(PBS_Seurat, file = ("PBS_Seurat_QC.RDS"))

R848_Seurat <- perform_qc(R848_Seurat)
saveRDS(R848_Seurat, file = ("R848_Seurat_QC.RDS"))

# Merge the Seurat objects
merged_seurat <- merge(Flag_Seurat, y = c(M229_Seurat, PBS_Seurat, R848_Seurat), add.cell.ids = c("Flag", "M229", "PBS", "R848"))
pbmc<-merged_seurat
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca") + NoLegend()
DimPlot(pbmc, reduction = "pca", group.by = "condition")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.6)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
pbmc <- RunUMAP(pbmc, dims = 1:20)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
# Assuming 'pbmc' is your Seurat object and the layers need to be joined
pbmc<-JoinLayers(pbmc)
saveRDS(pbmc,file="CLP_ToUse.RDS")
seurat <- pbmc

#Visualize UMAP labeled by vaccine group
p1 <- DimPlot(seurat, reduction = "umap", group.by = "condition")#+NoLegend()
p2 <- DimPlot(seurat, reduction = "umap", split.by = "condition")
plot(p1)
plot(p2)

p <- FeaturePlot(seurat, features = c("MS4A1","CD3E","LAMP3","CD14"), combine = FALSE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)


## Generating HeatMap to label cell populations ## 
## hallmark genes for immune cell populations ## 
### ALSO MINIMIZE NUMBER OF CELLS IN PLOT FOR EASE OF VIEWING ###

### ensure ident to plot is active
Idents(seurat) <- seurat@meta.data$seurat_clusters
### captures count of cells for the ident with the fewest
maxcells  <- min(table(Idents(seurat)))
### Markers to plot###
t.g <- c("CD3D","CD3G","CD3E","CD28","LCK","TRAT1","TRAC","SKAP1")
t4.g <- c("CD4")
t8.g <- c('CD8B','CD8A')
tr.g <- c("FOXP3","CTLA4","IL2RA")
dc.g <- c('CCL17','CCL22','CD1B','CD1E','IDO1','LAMP3')
b.g <- c('CD19','CD79A','CD79B','TCL1A','MS4A1','CD72')
pr.g <- c("MKI67",'CD70',"CD69","TYMS")
m.g <- c("CD14","CD163","S100A8","LYZ",'FCGR3A')
nk.g <- c('XCL1','XCL2','GZMB','KLRC1')
pdc.g <- c('IRF7', "IRF8", "TCF4", "FLT3", "CLEC4C", "LILRA4", 'IL3RA', "ZBTB46")
apop.g <- c('AIFM3','BAX','BCL10','CASP1','CASP4','BIK')
gly.g <- c('ENO1','GCK')
oxpho.g <- c('NDUFAF7','UQCRC1','SDHA','COA1','TRAP1','COX8A','COX7C','UQCR11','NDUFA13')
mit.g <- c('MTRNR2L10','MTRNR2L8','MTRNR2L3','MTRNR2L6','MTRNR2L11','IFNB2','HGF','CDF','FPRL1','CNTFR','WSX1')
mark1 <- c(t.g,t4.g,t8.g,tr.g,b.g,dc.g,m.g,nk.g,pr.g, apop.g, mit.g, pdc.g)
### nested object subsetting with downsampling
DoHeatmap(subset(seurat, downsample = maxcells), features = mark1, raster = FALSE)
ggsave('NoDual_CLP_MergedUMAP_CellTypeMarkers_Heatmap_byCluster.pdf',height=10,width=13, dpi = 300)

##Find top 5 Markers that distinguish each cluster from others##
seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE)

seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
DoHeatmap(subset(seurat, downsample = maxcells), features = top5$gene, size = 3, raster = FALSE) + NoLegend()
ggsave('NoDual_CLP_MergedUMAP_Top5Feature_Heatmap_byCluster.pdf',height=10,width=13, dpi = 300)

seurat$condition <- factor(x = seurat$condition, levels = c("PBS", "M229", "Flag", "R848"))
##Violin Plots with hallmark genes to visualize another way##
VlnPlot(seurat, group.by = "seurat_clusters", features = c(t.g, t4.g, t8.g, tr.g), pt.size = 0)
VlnPlot(seurat, group.by = "seurat_clusters", features = c(b.g), pt.size =0)
VlnPlot(seurat, group.by = "seurat_clusters", features = c(dc.g,m.g,nk.g), pt.size = 0)


#' -------- identified cell type labeling ----------
tem1 <- seurat@meta.data$seurat_clusters
annotation1 <- c('Naive CD4 T cells','ISG-stimulated cells','Naive CD8 T cells','B cells','Regulatory T cells','Effector CD4 T cells','Proliferating cells','Dendritic cells','Humanin+ cells','NK cells','Monocytes/Macs')

library(plyr)
for ( i in 0:10){
  tem2 <- mapvalues(tem1,i,annotation1[i+1]); tem1 <- tem2
}
seurat@meta.data$celltype <- as.character(tem2)

tem1 <- seurat@meta.data$seurat_clusters
order1 <- unique(as.character(na.omit(annotation1)));
col2 <- c("slategray2","violetred","palegreen3","royalblue3","goldenrod1", "orchid4","turquoise2","purple1","lightgrey","red1","sandybrown")

library(plyr); tem1 <- seurat@meta.data$celltype
for ( i in 1:length(order1)){
  tem2 <- mapvalues(tem1,order1[i],col2[i])
  tem1 <- tem2
}
seurat@meta.data$colors <- as.character(tem1)

library(ggplot2); library(RColorBrewer); library(ggsci)

order1 <- unique(as.character(na.omit(annotation1)));
seurat@meta.data$celltype <- factor(seurat@meta.data$celltype, levels=order1)
#clear environment of everything except combinedUMAP at this point - will take a lot of memory to run this code##
set.seed(1234)
DimPlot(seurat, reduction = "umap", group.by='celltype',label.size=2, cols=col2, pt.size = 0.3) 
rm('annotation1','col2','i','order1','samples','tem1','tem2')
saveRDS(seurat, file = "CLP_NoDual_mergedUMAPwCellTypes.RDS")