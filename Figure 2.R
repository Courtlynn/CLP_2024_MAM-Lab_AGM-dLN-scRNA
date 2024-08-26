
rm(list=ls())
gc()
#Load libraries

library(Seurat)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(Matrix)
library(dplyr)
library(readr)
library(glmGamPoi)

##set working directory where raw sequencing files are stored##
setwd("~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy")

##Load in Seurat Object with labeled cell types that was generated in figure 1##
##THIS IS MAIN OBJECT FOR REMAINDER OF STUDY##
seurat <- readRDS('CLP_NoDual_mergedUMAPwCellTypes.RDS')
Idents(seurat) <- seurat@meta.data$condition
samples <- c('PBS', 'M229', 'R848', 'Flag')
col2 <- c("slategray2","violetred","palegreen3","royalblue3","goldenrod1", "orchid4","turquoise2","purple1","lightgrey","red1","sandybrown")
##Splitting UMAP by Vaccine group - 4 UMAPs generated here##
set.seed(1234)
pst <- lapply(1:4,function(i){
  foo <- subset(seurat, subset = condition==samples[i])
  ps <- DimPlot(foo, reduction = "umap", label = FALSE, label.size=0, group.by = "celltype", cols=col2, pt.size = 0.3, alpha = 0.7)+ggtitle(samples[i])
  return (ps)})

library(cowplot)
pdf('NoDual_CLP_MergedUMAP_split_celltype.pdf',height=14,width=14)
plot_grid(plotlist=pst,ncol=2)
graphics.off()

#calculating number of cells per cluster per condition 
table(seurat@meta.data$celltype, seurat@meta.data$condition)
