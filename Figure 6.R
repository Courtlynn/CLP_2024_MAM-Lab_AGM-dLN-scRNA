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
#' ---------- create GO plot from sig. up/down gene lists ----------- '#
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
####### ISG-stimulated T cells ##########
######## IPR8+TLRa versus M229 ##########
#########################################

##############################
#### IPR8+flg v IPR8+m229 ###
##############################

ISG_FvM <- FindMarkers(seurat, ident.1 = "ISG-stimulated cells_Flag", ident.2 = "ISG-stimulated cells_M229", verbose = FALSE)
write.csv(ISG_FvM,file='ISG_FlagvM229.csv',quote=F)
## Remove Y linked genes ##
filtered_ISG_FvM_genes <- ISG_FvM[!rownames(ISG_FvM) %in% y_hgnc_symbol, ]
##Remove polymorphic HLA-II genes##
final_ISG_FvM_genes <- filtered_ISG_FvM_genes[!rownames(filtered_ISG_FvM_genes) %in% MHCII, ]
write.csv(final_ISG_FvM_genes,file='NoY_noMHCII_ISG_FlagvM229.csv',quote=F)

######################################################################
## Build the transcript ID to gene ID table from the entrez database##
######################################################################

humanannot<-grch38
humanannot
dfg <- final_ISG_FvM_genes
dfg <- as.data.frame(dfg)
dfg$RowNames <- rownames(dfg)
colnames(dfg)[6] <- "RowNames"
df <- cbind(humanannot$entrez, humanannot$symbol)
df <- as.data.frame(df)
dfg
df
final_ISG_FvM_1 <- left_join(dfg,df,by=c("RowNames"="V2"))
View(final_ISG_FvM_1)

########################################################
### Generate list of significantly UPERGULATED genes ###
########################################################

final_ISG_FvM_1 <- final_ISG_FvM_1[final_ISG_FvM_1$p_val_adj < 0.05,]
final_ISG_FvM_2 <- final_ISG_FvM_1[final_ISG_FvM_1$avg_log2FC > 0,]
gene_list_up = final_ISG_FvM_2$avg_log2FC
names(gene_list_up) = final_ISG_FvM_2$V1
gene_list_up = sort(gene_list_up, decreasing = TRUE)
gene_list_up = gene_list_up[!duplicated(names(gene_list_up))] #removes duplicates based on name
#gene_list = gene_list[!duplicated(gene_list)] #removes duplicates based on value - should we include this???
head(gene_list_up,n=500)

##ClusterProfiler GO enrichment##
##UPREGULATED GENE LIST###
set.seed(1234)
GO_ISG_FvM_up<- enrichGO(names(gene_list_up), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_ISG_FvM_up_1 <- clusterProfiler:: simplify(GO_ISG_FvM_up, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_ISG_FvM_up_1 <- dotplot(GO_ISG_FvM_up_1)+ ggtitle(label = "Up_GO_Enrichment: IPR8+flg v IPR8+m229", subtitle = "ISG-stimulated T cells")
dotplt_ISG_FvM_up_1 ##Export plot as PNG

##########################################################
### Generate list of significantly DOWNREGULATED genes ###
##########################################################

final_ISG_FvM_1 <- final_ISG_FvM_1[final_ISG_FvM_1$p_val_adj < 0.05,]
final_ISG_FvM_3 <- final_ISG_FvM_1[final_ISG_FvM_1$avg_log2FC < 0,]
gene_list_dwn = final_ISG_FvM_3$avg_log2FC
names(gene_list_dwn) = final_ISG_FvM_3$V1
gene_list_dwn = sort(gene_list_dwn, decreasing = TRUE)
gene_list_dwn = gene_list_dwn[!duplicated(names(gene_list_dwn))] #removes duplicates based on name
#gene_list = gene_list[!duplicated(gene_list)] #removes duplicates based on value - should we include this???
head(gene_list_dwn,n=500)

###DOWNREGULATED GENE LIST###
set.seed(1234)
GO_ISG_FvM_dwn<- enrichGO(names(gene_list_dwn), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_ISG_FvM_dwn_1 <- clusterProfiler:: simplify(GO_ISG_FvM_dwn, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_ISG_FvM_dwn_1 <- dotplot(GO_ISG_FvM_dwn_1)+ ggtitle(label = "Down_GO_Enrichment: IPR8+flg v IPR8+m229", subtitle = "ISG-stimulated T cells")
dotplt_ISG_FvM_dwn_1 ##Export plot as PNG

##############################
#### IPR8-R848 v IPR8+m229 ###
##############################

ISG_RvM <- FindMarkers(seurat, ident.1 = "ISG-stimulated cells_R848", ident.2 = "ISG-stimulated cells_M229", verbose = FALSE)
write.csv(ISG_RvM,file='ISG_R848vM229.csv',quote=F)
## Remove Y linked genes ##
filtered_ISG_RvM_genes <- ISG_RvM[!rownames(ISG_RvM) %in% y_hgnc_symbol, ]
##Remove polymorphic HLA-II genes##
final_ISG_RvM_genes <- filtered_ISG_RvM_genes[!rownames(filtered_ISG_RvM_genes) %in% MHCII, ]
write.csv(final_ISG_RvM_genes,file='NoY_noMHCII_ISG_R848vM229.csv',quote=F)

######################################################################
## Build the transcript ID to gene ID table from the entrez database##
######################################################################

humanannot<-grch38
humanannot
dfg <- final_ISG_RvM_genes
dfg <- as.data.frame(dfg)
dfg$RowNames <- rownames(dfg)
colnames(dfg)[6] <- "RowNames"
df <- cbind(humanannot$entrez, humanannot$symbol)
df <- as.data.frame(df)
dfg
df
final_ISG_RvM_1 <- left_join(dfg,df,by=c("RowNames"="V2"))
View(final_ISG_RvM_1)

########################################################
### Generate list of significantly UPERGULATED genes ###
########################################################

final_ISG_RvM_1 <- final_ISG_RvM_1[final_ISG_RvM_1$p_val_adj < 0.05,]
final_ISG_RvM_2 <- final_ISG_RvM_1[final_ISG_RvM_1$avg_log2FC > 0,]
gene_list_up = final_ISG_RvM_2$avg_log2FC
names(gene_list_up) = final_ISG_RvM_2$V1
gene_list_up = sort(gene_list_up, decreasing = TRUE)
gene_list_up = gene_list_up[!duplicated(names(gene_list_up))] #removes duplicates based on name
#gene_list = gene_list[!duplicated(gene_list)] #removes duplicates based on value - should we include this???
head(gene_list_up,n=500)

##ClusterProfiler GO enrichment##
##UPREGULATED GENE LIST###
set.seed(1234)
GO_ISG_RvM_up<- enrichGO(names(gene_list_up), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_ISG_RvM_up_1 <- clusterProfiler:: simplify(GO_ISG_RvM_up, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_ISG_RvM_up_1 <- dotplot(GO_ISG_RvM_up_1)+ ggtitle(label = "Up_GO_Enrichment: IPR8-R848 v IPR8+m229", subtitle = "ISG-stimulated T cells")
dotplt_ISG_RvM_up_1 ##Export plot as PNG

##########################################################
### Generate list of significantly DOWNREGULATED genes ###
##########################################################

final_ISG_RvM_1 <- final_ISG_RvM_1[final_ISG_RvM_1$p_val_adj < 0.05,]
final_ISG_RvM_3 <- final_ISG_RvM_1[final_ISG_RvM_1$avg_log2FC < 0,]
gene_list_dwn = final_ISG_RvM_3$avg_log2FC
names(gene_list_dwn) = final_ISG_RvM_3$V1
gene_list_dwn = sort(gene_list_dwn, decreasing = TRUE)
gene_list_dwn = gene_list_dwn[!duplicated(names(gene_list_dwn))] #removes duplicates based on name
#gene_list = gene_list[!duplicated(gene_list)] #removes duplicates based on value - should we include this???
head(gene_list_dwn,n=500)

###DOWNREGULATED GENE LIST###
set.seed(1234)
GO_ISG_RvM_dwn<- enrichGO(names(gene_list_dwn), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_ISG_RvM_dwn_1 <- clusterProfiler:: simplify(GO_ISG_RvM_dwn, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_ISG_RvM_dwn_1 <- dotplot(GO_ISG_RvM_dwn_1)+ ggtitle(label = "Down_GO_Enrichment: IPR8-R848 v IPR8+m229", subtitle = "ISG-stimulated T cells")
dotplt_ISG_RvM_dwn_1 ##Export plot as PNG

##############################
#### IPR8-R848 v IPR8+flg ####
##############################

ISG_RvF <- FindMarkers(seurat, ident.1 = "ISG-stimulated cells_R848", ident.2 = "ISG-stimulated cells_Flag", verbose = FALSE)
write.csv(ISG_RvF,file='ISG_R848vFlag.csv',quote=F)
## Remove Y linked genes ##
filtered_ISG_RvF_genes <- ISG_RvF[!rownames(ISG_RvF) %in% y_hgnc_symbol, ]
##Remove polymorphic HLA-II genes##
final_ISG_RvF_genes <- filtered_ISG_RvF_genes[!rownames(filtered_ISG_RvF_genes) %in% MHCII, ]
write.csv(final_ISG_RvF_genes,file='NoY_noMHCII_ISG_R848vFlag.csv',quote=F)

######################################################################
## Build the transcript ID to gene ID table from the entrez database##
######################################################################

humanannot<-grch38
humanannot
dfg <- final_ISG_RvF_genes
dfg <- as.data.frame(dfg)
dfg$RowNames <- rownames(dfg)
colnames(dfg)[6] <- "RowNames"
df <- cbind(humanannot$entrez, humanannot$symbol)
df <- as.data.frame(df)
dfg
df
final_ISG_RvF_1 <- left_join(dfg,df,by=c("RowNames"="V2"))
View(final_ISG_RvF_1)

########################################################
### Generate list of significantly UPERGULATED genes ###
########################################################

final_ISG_RvF_1 <- final_ISG_RvF_1[final_ISG_RvF_1$p_val_adj < 0.05,]
final_ISG_RvF_2 <- final_ISG_RvF_1[final_ISG_RvF_1$avg_log2FC > 0,]
gene_list_up = final_ISG_RvF_2$avg_log2FC
names(gene_list_up) = final_ISG_RvF_2$V1
gene_list_up = sort(gene_list_up, decreasing = TRUE)
gene_list_up = gene_list_up[!duplicated(names(gene_list_up))] #removes duplicates based on name
#gene_list = gene_list[!duplicated(gene_list)] #removes duplicates based on value - should we include this???
head(gene_list_up,n=500)

##ClusterProfiler GO enrichment##
##UPREGULATED GENE LIST###
set.seed(1234)
GO_ISG_RvF_up<- enrichGO(names(gene_list_up), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_ISG_RvF_up_1 <- clusterProfiler:: simplify(GO_ISG_RvF_up, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_ISG_RvF_up_1 <- dotplot(GO_ISG_RvF_up_1)+ ggtitle(label = "Up_GO_Enrichment: IPR8-R848 v IPR8+flg", subtitle = "ISG-stimulated T cells")
dotplt_ISG_RvF_up_1 ##Export plot as PNG

##########################################################
### Generate list of significantly DOWNREGULATED genes ###
##########################################################

final_ISG_RvF_1 <- final_ISG_RvF_1[final_ISG_RvF_1$p_val_adj < 0.05,]
final_ISG_RvF_3 <- final_ISG_RvF_1[final_ISG_RvF_1$avg_log2FC < 0,]
gene_list_dwn = final_ISG_RvF_3$avg_log2FC
names(gene_list_dwn) = final_ISG_RvF_3$V1
gene_list_dwn = sort(gene_list_dwn, decreasing = TRUE)
gene_list_dwn = gene_list_dwn[!duplicated(names(gene_list_dwn))] #removes duplicates based on name
#gene_list = gene_list[!duplicated(gene_list)] #removes duplicates based on value - should we include this???
head(gene_list_dwn,n=500)

###DOWNREGULATED GENE LIST###
set.seed(1234)
GO_ISG_RvF_dwn<- enrichGO(names(gene_list_dwn), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_ISG_RvF_dwn_1 <- clusterProfiler:: simplify(GO_ISG_RvF_dwn, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_ISG_RvF_dwn_1 <- dotplot(GO_ISG_RvF_dwn_1)+ ggtitle(label = "Down_GO_Enrichment: IPR8-R848 v IPR8+flg", subtitle = "ISG-stimulated T cells")
dotplt_ISG_RvF_dwn_1 ##Export plot as PNG

##Top 25 most signficantly upregulated genes##
ISGtop25_RvF <- final_ISG_RvF_genes %>%  
  dplyr::filter(p_val_adj < 0.05) %>% #only significant genes
  dplyr::filter(avg_log2FC > 0) %>% #take upregulated genes
  dplyr::arrange(p_val_adj) %>%  # Sort by p_adj_val 
  slice_head(n = 25) %>%
  ungroup() -> ISGtop25_RvF

Idents(seurat) <- seurat@meta.data$celltype
DotPlot(seurat, features = rownames(ISGtop25_RvF), scale = T, dot.scale = 10, group.by = "condition" ,idents = "ISG-stimulated cells")+
  theme_minimal() + ggtitle(label = "ISG-stimulated T cells: IPR8-R848 v IPR8+flg", subtitle = "Top 25 most signifcant DEG (Log2FC > 0)") + # Change the theme to minimal and add title 
  theme(plot.background = element_rect(fill = "white")) +  # Set plot background to white
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, margin = margin(t = 10))) +  # Rotate x-axis labels by 45 degrees and adjust horizontal justification
  scale_colour_gradient2(low="black", mid="lightgrey", high="darkgoldenrod")

##Top 25 most significantly downregulated genes##
ISGbot25_RvF <- final_ISG_RvF_genes %>%  
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr:: filter(avg_log2FC < -0) %>%
  dplyr:: arrange(p_val_adj)%>%
  slice_head(n = 25) %>%
  ungroup() -> ISGbot25_RvF

Idents(seurat) <- seurat@meta.data$celltype
DotPlot(seurat, features = rownames(ISGbot25_RvF), scale = T, dot.scale = 10, group.by = "condition", idents = "ISG-stimulated cells")+
  theme_minimal() + ggtitle(label = "ISG-stimulated T cells: IPR8-R848 v IPR8+flg", subtitle = "Top 25 most signifcant DEG (Log2FC < 0)") +  # Change the theme to minimal
  theme(plot.background = element_rect(fill = "white")) +  # Set plot background to white
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, margin = margin(t = 10))) +  # Rotate x-axis labels by 45 degrees and adjust horizontal justification
  scale_colour_gradient2(low="black", mid="lightgrey", high="darkgoldenrod")
