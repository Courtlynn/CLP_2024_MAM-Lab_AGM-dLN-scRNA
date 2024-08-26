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
library(biomaRt)
library(data.table)
library(cowplot)
library(ggplot2)
library(tidyverse)
library(devtools)
library(clusterProfiler)
library(mygene) 

##set working directory where raw sequencing files are stored##
setwd("~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy")

##Load in Seurat Object with labeled cell types that was generated in figure 1##
##THIS IS MAIN OBJECT FOR REMAINDER OF STUDY##
seurat <- readRDS('CLP_NoDual_mergedUMAPwCellTypes.RDS')

##############################################################
############### Overview of following code ###################
##############################################################

#' ----- DEGs for IPR8+TLRa v IPR8+m229 and IPR8-R848 v IPR8+flg ---- '#
#' ---------- remove Y-linked and polymorphic genes ----------- '#
#' ------ create list of signifcantly up and downregulated genes ----- '#
#' ------- create Dotplot of top 25 DEGs -------- '#

## Save DEG analysis as csv file after each run ##

##################################################################
### Y linked genes will need to be removed from DEG lists###
########### First make list of all Y linked genes###############
##################################################################

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
y_genes <- getBM(attributes = c("chromosome_name", "ensembl_gene_id", "hgnc_symbol"),
                 filters = "chromosome_name", values = c("Y"), mart = mart)# count genes

y_hgnc_symbol <- y_genes$hgnc_symbol

##MHCII genes are polymorphic and therefore will need to be removed after DEG generation##
MHCII <- c("HLA-DRB5", "HLA-DRB1", "HLA-DRA", "HLA-DPB1", "HLA-DPA1",  "HLA-DQA1", "HLA-DQA2","HLA-DQB2")

###################################
### APCs: IPR8+flg v IPR8+m229 ####
###################################
APC_FvM <- FindMarkers(seurat, ident.1 = c("Dendritic cells_Flag", "Monocytes/Macs_Flag"), ident.2 = c("Dendritic cells_M229", "Monocytes/Macs_M229"), verbose = FALSE)
write.csv(APC_FvM,file='APCs_FlagvM229.csv',quote=F)
## Remove Y linked genes ##
filtered_APC_FvM_genes <- APC_FvM[!rownames(APC_FvM) %in% y_hgnc_symbol, ]
##Remove polymorphic HLA-II genes##
final_APC_FvM_genes <- filtered_APC_FvM_genes[!rownames(filtered_APC_FvM_genes) %in% MHCII, ]
write.csv(final_APC_FvM_genes,file='NoY_noMHCII_APC_FlagvM229.csv',quote=F)

##Top 25 most signficantly upregulated genes##
APCtop25_FvM <- final_APC_FvM_genes %>%  
  dplyr::filter(p_val_adj < 0.05) %>% #only significant genes
  dplyr::filter(avg_log2FC > 0) %>% #take upregulated genes
  dplyr::arrange(p_val_adj) %>%  # Sort by p_adj_val 
  slice_head(n = 25) %>%
  ungroup() -> APCtop25_FvM

Idents(seurat) <- seurat@meta.data$celltype
DotPlot(seurat, features = rownames(APCtop25_FvM), scale = T, dot.scale = 10, group.by = "condition", idents = c("Dendritic cells", "Monocytes/Macs"))+
  theme_minimal() + ggtitle(label = "APCs: IPR8+flg vs IPR8+m229", subtitle = "Top 25 most signifcant DEG (Log2FC > 0)") + # Change the theme to minimal and add title 
  theme(plot.background = element_rect(fill = "white")) +  # Set plot background to white
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, margin = margin(t = 10))) +  # Rotate x-axis labels by 45 degrees and adjust horizontal justification
  scale_colour_gradient2(low="black", mid="lightgrey", high="darkgoldenrod")

##Top 25 most significantly downregulated genes##
APCbot25_FvM <- final_APC_FvM_genes %>%  
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr:: filter(avg_log2FC < -0) %>%
  dplyr:: arrange(p_val_adj)%>%
  slice_head(n = 25) %>%
  ungroup() -> APCbot25_FvM

Idents(seurat) <- seurat@meta.data$celltype
DotPlot(seurat, features = rownames(APCbot25_FvM), scale = T, dot.scale = 10, group.by = "condition", idents = c("Dendritic cells", "Monocytes/Macs"))+
  theme_minimal() + ggtitle(label = "APCs: IPR8+flg vs IPR8+m229", subtitle = "Top 25 most signifcant DEG (Log2FC < 0)") +  # Change the theme to minimal
  theme(plot.background = element_rect(fill = "white")) +  # Set plot background to white
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, margin = margin(t = 10))) +  # Rotate x-axis labels by 45 degrees and adjust horizontal justification
  scale_colour_gradient2(low="black", mid="lightgrey", high="darkgoldenrod")


###################################
### APCs: IPR8-R848 v IPR8+m229 ####
###################################
APC_RvM <- FindMarkers(seurat, ident.1 = c("Dendritic cells_R848", "Monocytes/Macs_R848"), ident.2 = c("Dendritic cells_M229", "Monocytes/Macs_M229"), verbose = FALSE)
write.csv(APC_RvM,file='APCs_R848vM229.csv',quote=F)
## Remove Y linked genes ##
filtered_APC_RvM_genes <- APC_RvM[!rownames(APC_RvM) %in% y_hgnc_symbol, ]
##Remove polymorphic HLA-II genes##
final_APC_RvM_genes <- filtered_APC_RvM_genes[!rownames(filtered_APC_RvM_genes) %in% MHCII, ]
write.csv(final_APC_RvM_genes,file='NoY_noMHCII_APC_R848vM229.csv',quote=F)

##Top 25 most signficantly upregulated genes##
APCtop25_RvM <- final_APC_RvM_genes %>%  
  dplyr::filter(p_val_adj < 0.05) %>% #only significant genes
  dplyr::filter(avg_log2FC > 0) %>% #take upregulated genes
  dplyr::arrange(p_val_adj) %>%  # Sort by p_adj_val 
  slice_head(n = 25) %>%
  ungroup() -> APCtop25_RvM

Idents(seurat) <- seurat@meta.data$celltype
DotPlot(seurat, features = rownames(APCtop25_RvM), scale = T, dot.scale = 10, group.by = "condition", idents = c("Dendritic cells", "Monocytes/Macs"))+
  theme_minimal() + ggtitle(label = "APCs: IPR8-R848 vs IPR8+m229", subtitle = "Top 25 most signifcant DEG (Log2FC > 0)") + # Change the theme to minimal and add title 
  theme(plot.background = element_rect(fill = "white")) +  # Set plot background to white
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, margin = margin(t = 10))) +  # Rotate x-axis labels by 45 degrees and adjust horizontal justification
  scale_colour_gradient2(low="black", mid="lightgrey", high="darkgoldenrod")

##Top 25 most significantly downregulated genes##
APCbot25_RvM <- final_APC_RvM_genes %>%  
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr:: filter(avg_log2FC < -0) %>%
  dplyr:: arrange(p_val_adj)%>%
  slice_head(n = 25) %>%
  ungroup() -> APCbot25_RvM

Idents(seurat) <- seurat@meta.data$celltype
DotPlot(seurat, features = rownames(APCbot25_RvM), scale = T, dot.scale = 10, group.by = "condition", idents = c("Dendritic cells", "Monocytes/Macs"))+
  theme_minimal() + ggtitle(label = "APCs: IPR8-R848 vs IPR8+m229", subtitle = "Top 25 most signifcant DEG (Log2FC < 0)") +  # Change the theme to minimal
  theme(plot.background = element_rect(fill = "white")) +  # Set plot background to white
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, margin = margin(t = 10))) +  # Rotate x-axis labels by 45 degrees and adjust horizontal justification
  scale_colour_gradient2(low="black", mid="lightgrey", high="darkgoldenrod")

###################################
### APCs: IPR8-R848 v IPR8+flg ####
###################################
APC_RvF <- FindMarkers(seurat, ident.1 = c("Dendritic cells_R848", "Monocytes/Macs_R848"), ident.2 = c("Dendritic cells_Flag", "Monocytes/Macs_Flag"), verbose = FALSE)
write.csv(APC_RvF,file='APCs_R848vFlag.csv',quote=F)
## Remove Y linked genes ##
filtered_APC_RvF_genes <- APC_RvF[!rownames(APC_RvF) %in% y_hgnc_symbol, ]
##Remove polymorphic HLA-II genes##
final_APC_RvF_genes <- filtered_APC_RvF_genes[!rownames(filtered_APC_RvF_genes) %in% MHCII, ]
write.csv(final_APC_RvF_genes,file='NoY_noMHCII_APC_R848vFlag.csv',quote=F)


##Top 25 most signficantly upregulated genes##
APCtop25_RvF <- final_APC_RvF_genes %>%  
  dplyr::filter(p_val_adj < 0.05) %>% #only significant genes
  dplyr::filter(avg_log2FC > 0) %>% #take upregulated genes
  dplyr::arrange(p_val_adj) %>%  # Sort by p_adj_val 
  slice_head(n = 25) %>%
  ungroup() -> APCtop25_RvF

Idents(seurat) <- seurat@meta.data$celltype
DotPlot(seurat, features = rownames(APCtop25_RvF), scale = T, dot.scale = 10, group.by = "condition", idents = c("Dendritic cells", "Monocytes/Macs"))+
  theme_minimal() + ggtitle(label = "APCs: IPR8-R848 vs IPR8+flg", subtitle = "Top 25 most signifcant DEG (Log2FC > 0)") + # Change the theme to minimal and add title 
  theme(plot.background = element_rect(fill = "white")) +  # Set plot background to white
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, margin = margin(t = 10))) +  # Rotate x-axis labels by 45 degrees and adjust horizontal justification
  scale_colour_gradient2(low="black", mid="lightgrey", high="darkgoldenrod")

##Top 25 most significantly downregulated genes##
APCbot25_RvF <- final_APC_RvF_genes %>%  
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr:: filter(avg_log2FC < -0) %>%
  dplyr:: arrange(p_val_adj)%>%
  slice_head(n = 25) %>%
  ungroup() -> APCbot25_RvF

Idents(seurat) <- seurat@meta.data$celltype
DotPlot(seurat, features = rownames(APCbot25_RvF), scale = T, dot.scale = 10, group.by = "condition", idents = c("Dendritic cells", "Monocytes/Macs"))+
  theme_minimal() + ggtitle(label = "APCs: IPR8-R848 vs IPR8+flg", subtitle = "Top 25 most signifcant DEG (Log2FC < 0)") +  # Change the theme to minimal
  theme(plot.background = element_rect(fill = "white")) +  # Set plot background to white
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, margin = margin(t = 10))) +  # Rotate x-axis labels by 45 degrees and adjust horizontal justification
  scale_colour_gradient2(low="black", mid="lightgrey", high="darkgoldenrod")
