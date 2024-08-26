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

#' ----- DEGs for IPR8+TLRa v IPR8+m229, then IPR8-R848 v IPR8+flg ---- '#
#' ---------- remove Y-linked and polymorphic genes ----------- '#
#' ------ create list of signifcantly up and downregulated genes ----- '#
#' ------- create Dotplot of Top25 sig up/down gene lists (RvF) -------- '#

## For the All Clusters analyses below, min.pct is set to 0.01 ##
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

#########################################
############# All Cells #################
######## IPR8+TLRa versus M229 ##########
#########################################

##############################
#### IPR8-R848 v IPR8+m229 ###
##############################

RvM <- FindMarkers(seurat,  ident.1 = 'R848', ident.2 = 'M229', only.pos = FALSE, min.pct = 0.01, logfc.threshold = -Inf, min.cells.feature = 1, min.cells.group = 1)
write.csv(RvM,file='01_AllClusters_R848vM229.csv',quote=F)
## Remove Y linked genes ##
filtered_RvM_genes <- RvM[!rownames(RvM) %in% y_hgnc_symbol, ]
##Remove polymorphic HLA-II genes##
final_RvM_genes <- filtered_RvM_genes[!rownames(filtered_RvM_genes) %in% MHCII, ]
write.csv(final_RvM_genes,file='NoY_noMHCII_01_AllClusters_R848vM229.csv',quote=F)

##Top 25 most signficantly upregulated genes##
Alltop25_RvM <- final_RvM_genes %>%  
  dplyr::filter(p_val_adj < 0.05) %>% #only significant genes
  dplyr::filter(avg_log2FC > 0) %>% #take upregulated genes
  dplyr::arrange(p_val_adj) %>%  # Sort by p_adj_val 
  slice_head(n = 25) %>%
  ungroup() -> Alltop25_RvM

DotPlot(seurat, features = rownames(Alltop25_RvM), scale = T, dot.scale = 10, group.by = "condition")+
  theme_minimal() + ggtitle(label = "All cells: IPR8-R848 v IPR8+m229", subtitle = "Top 25 most signifcant DEG (Log2FC > 0)") + # Change the theme to minimal and add title 
  theme(plot.background = element_rect(fill = "white")) +  # Set plot background to white
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, margin = margin(t = 10))) +  # Rotate x-axis labels by 45 degrees and adjust horizontal justification
  scale_colour_gradient2(low="black", mid="lightgrey", high="darkgoldenrod")


##Top 25 most significantly downregulated genes##
Allbot25_RvM <- final_RvM_genes %>%  
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr:: filter(avg_log2FC < -0) %>%
  dplyr:: arrange(p_val_adj)%>%
  slice_head(n = 25) %>%
  ungroup() -> Allbot25_RvM


DotPlot(seurat, features = rownames(Allbot25_RvM), scale = T, dot.scale = 10, group.by = "condition")+
  theme_minimal() + ggtitle(label = "All cells: IPR8-R848 v IPR8+m229", subtitle = "Top 25 most signifcant DEG (Log2FC < 0)") +  # Change the theme to minimal
  theme(plot.background = element_rect(fill = "white")) +  # Set plot background to white
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, margin = margin(t = 10))) +  # Rotate x-axis labels by 45 degrees and adjust horizontal justification
  scale_colour_gradient2(low="black", mid="lightgrey", high="darkgoldenrod")

##############################
#### IPR8+flg v IPR8+m229 ###
##############################

FvM <- FindMarkers(seurat, ident.1 = 'Flag', ident.2 = 'M229', only.pos = FALSE, min.pct = 0.01, logfc.threshold = -Inf, min.cells.feature = 1, min.cells.group = 1)
write.csv(FvM,file='01_AllClusters_FlagvM229.csv',quote=F)
## Remove X and Y linked genes ##
filtered_FvM_genes <- FvM[!rownames(FvM) %in% y_hgnc_symbol, ]
##Remove polymorphic HLA-II genes##
final_FvM_genes <- filtered_FvM_genes[!rownames(filtered_FvM_genes) %in% MHCII, ]
write.csv(final_FvM_genes,file='NoY_noMHCII_01_AllClusters_FlagvM229.csv',quote=F)

##Top 25 most signficantly upregulated genes##
Alltop25_FvM <- final_FvM_genes %>%  
  dplyr::filter(p_val_adj < 0.05) %>% #only significant genes
  dplyr::filter(avg_log2FC > 0) %>% #take upregulated genes
  dplyr::arrange(p_val_adj) %>%  # Sort by p_adj_val 
  slice_head(n = 25) %>%
  ungroup() -> Alltop25_FvM

DotPlot(seurat, features = rownames(Alltop25_FvM), scale = T, dot.scale = 10, group.by = "condition")+
  theme_minimal() + ggtitle(label = "All cells: IPR8+flg v IPR8+m229", subtitle = "Top 25 most signifcant DEG (Log2FC > 0)") + # Change the theme to minimal and add title 
  theme(plot.background = element_rect(fill = "white")) +  # Set plot background to white
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, margin = margin(t = 10))) +  # Rotate x-axis labels by 45 degrees and adjust horizontal justification
  scale_colour_gradient2(low="black", mid="lightgrey", high="darkgoldenrod")


##Top 25 most significantly downregulated genes##
Allbot25_FvM <- final_FvM_genes %>%  
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr:: filter(avg_log2FC < -0) %>%
  dplyr:: arrange(p_val_adj)%>%
  slice_head(n = 25) %>%
  ungroup() -> Allbot25_FvM


DotPlot(seurat, features = rownames(Allbot25_FvM), scale = T, dot.scale = 10, group.by = "condition")+
  theme_minimal() + ggtitle(label = "All cells: IPR8+flg v IPR8+m229", subtitle = "Top 25 most signifcant DEG (Log2FC < 0)") +  # Change the theme to minimal
  theme(plot.background = element_rect(fill = "white")) +  # Set plot background to white
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, margin = margin(t = 10))) +  # Rotate x-axis labels by 45 degrees and adjust horizontal justification
  scale_colour_gradient2(low="black", mid="lightgrey", high="darkgoldenrod")