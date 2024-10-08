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


######################################################################
## Build the transcript ID to gene ID table from the entrez database##
######################################################################

humanannot<-grch38
humanannot
dfg <- final_RvM_genes
dfg <- as.data.frame(dfg)
dfg$RowNames <- rownames(dfg)
colnames(dfg)[6] <- "RowNames"
df <- cbind(humanannot$entrez, humanannot$symbol)
df <- as.data.frame(df)
dfg
df
final_RvM_1 <- left_join(dfg,df,by=c("RowNames"="V2"))
View(final_RvM_1)

########################################################
### Generate list of significantly UPERGULATED genes ###
########################################################

final_RvM_1 <- final_RvM_1[final_RvM_1$p_val_adj < 0.05,]
final_RvM_2 <- final_RvM_1[final_RvM_1$avg_log2FC > 0,]
gene_list_up = final_RvM_2$avg_log2FC
names(gene_list_up) = final_RvM_2$V1
gene_list_up = sort(gene_list_up, decreasing = TRUE)
gene_list_up = gene_list_up[!duplicated(names(gene_list_up))] #removes duplicates based on name
head(gene_list_up,n=500)

##ClusterProfiler GO enrichment##
##UPREGULATED GENE LIST###
set.seed(1234)
GO_RvM_up<- enrichGO(names(gene_list_up), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_RvM_up_1 <- clusterProfiler:: simplify(GO_RvM_up, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_RvM_up_1 <- dotplot(GO_RvM_up_1)+ ggtitle(label = "Up_GO_Enrichment: IPR8-R848 v IPR8+m229", subtitle = "All Cells")
dotplt_RvM_up_1 ##Export Plot as PNG 

##########################################################
### Generate list of significantly DOWNREGULATED genes ###
##########################################################

final_RvM_1 <- final_RvM_1[final_RvM_1$p_val_adj < 0.05,]
final_RvM_3 <- final_RvM_1[final_RvM_1$avg_log2FC < 0,]
gene_list_dwn = final_RvM_3$avg_log2FC
names(gene_list_dwn) = final_RvM_3$V1
gene_list_dwn = sort(gene_list_dwn, decreasing = TRUE)
gene_list_dwn = gene_list_dwn[!duplicated(names(gene_list_dwn))] #removes duplicates based on name
head(gene_list_dwn,n=500)

###DOWNREGULATED GENE LIST###
set.seed(1234)
GO_RvM_dwn<- enrichGO(names(gene_list_dwn), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_RvM_dwn_1 <- clusterProfiler:: simplify(GO_RvM_dwn, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_RvM_dwn_1 <- dotplot(GO_RvM_dwn_1)+ ggtitle(label = "Down_GO_Enrichment: IPR8-R848 v IPR8+m229", subtitle = "All Cells")
dotplt_RvM_dwn_1 ##Export Plot as PNG 

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

######################################################################
## Build the transcript ID to gene ID table from the entrez database##
######################################################################

humanannot<-grch38
humanannot
dfg <- final_FvM_genes
dfg <- as.data.frame(dfg)
dfg$RowNames <- rownames(dfg)
colnames(dfg)[6] <- "RowNames"
df <- cbind(humanannot$entrez, humanannot$symbol)
df <- as.data.frame(df)
dfg
df
final_FvM_1 <- left_join(dfg,df,by=c("RowNames"="V2"))
View(final_FvM_1)

########################################################
### Generate list of significantly UPERGULATED genes ###
########################################################

final_FvM_1 <- final_FvM_1[final_FvM_1$p_val_adj < 0.05,]
final_FvM_2 <- final_FvM_1[final_FvM_1$avg_log2FC > 0,]
gene_list_up = final_FvM_2$avg_log2FC
names(gene_list_up) = final_FvM_2$V1
gene_list_up = sort(gene_list_up, decreasing = TRUE)
gene_list_up = gene_list_up[!duplicated(names(gene_list_up))] #removes duplicates based on name
head(gene_list_up,n=500)

##ClusterProfiler GO enrichment##
##UPREGULATED GENE LIST###
set.seed(1234)
GO_FvM_up<- enrichGO(names(gene_list_up), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_FvM_up_1 <- clusterProfiler:: simplify(GO_FvM_up, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_FvM_up_1 <- dotplot(GO_FvM_up_1)+ ggtitle(label = "Up_GO_Enrichment: IPR8+flg v IPR8+m229", subtitle = "All Cells")
dotplt_FvM_up_1 ##Export Plot as PNG 

##########################################################
### Generate list of significantly DOWNREGULATED genes ###
##########################################################

final_FvM_1 <- final_FvM_1[final_FvM_1$p_val_adj < 0.05,]
final_FvM_3 <- final_FvM_1[final_FvM_1$avg_log2FC < 0,]
gene_list_dwn = final_FvM_3$avg_log2FC
names(gene_list_dwn) = final_FvM_3$V1
gene_list_dwn = sort(gene_list_dwn, decreasing = TRUE)
gene_list_dwn = gene_list_dwn[!duplicated(names(gene_list_dwn))] #removes duplicates based on name
head(gene_list_dwn,n=500)

###DOWNREGULATED GENE LIST###
set.seed(1234)
GO_FvM_dwn<- enrichGO(names(gene_list_dwn), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_FvM_dwn_1 <- clusterProfiler:: simplify(GO_FvM_dwn, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_FvM_dwn_1 <- dotplot(GO_FvM_dwn_1)+ ggtitle(label = "Down_GO_Enrichment: IPR8+flg v IPR8+m229", subtitle = "All Cells")
dotplt_FvM_dwn_1 ##Export Plot as PNG 

#########################################
############# All Cells #################
######## IPR8-R848 v IPR8+flg ##########
#########################################

###################################################
## Finally, directly compare two adjuvant groups ##
###################################################
RvF <- FindMarkers(seurat, ident.1 = 'R848', ident.2 = 'Flag', only.pos = FALSE, min.pct = -Inf, logfc.threshold = -Inf, min.cells.feature = 1, min.cells.group = 1)
write.csv(RvF,file='01_AllClusters_R848vFlag.csv',quote=F)
## Remove X and Y linked genes ##
filtered_RvF_genes <- RvF[!rownames(RvF) %in% y_hgnc_symbol, ]
##Remove polymorphic HLA-II genes##
final_RvF_genes <- filtered_RvF_genes[!rownames(filtered_RvF_genes) %in% MHCII, ]
write.csv(final_RvF_genes,file='NoY_noMHCII_01_AllClusters_R848vFlag.csv',quote=F)

######################################################################
## Build the transcript ID to gene ID table from the entrez database##
######################################################################

humanannot<-grch38
humanannot
dfg <- final_RvF_genes
dfg <- as.data.frame(dfg)
dfg$RowNames <- rownames(dfg)
colnames(dfg)[6] <- "RowNames"
df <- cbind(humanannot$entrez, humanannot$symbol)
df <- as.data.frame(df)
dfg
df
final_RvF_1 <- left_join(dfg,df,by=c("RowNames"="V2"))
View(final_FvM_1)

########################################################
### Generate list of significantly UPERGULATED genes ###
########################################################

final_RvF_1 <- final_RvF_1[final_RvF_1$p_val_adj < 0.05,]
final_RvF_2 <- final_RvF_1[final_RvF_1$avg_log2FC > 0,]
gene_list_up = final_RvF_2$avg_log2FC
names(gene_list_up) = final_RvF_2$V1
gene_list_up = sort(gene_list_up, decreasing = TRUE)
gene_list_up = gene_list_up[!duplicated(names(gene_list_up))] #removes duplicates based on name
head(gene_list_up,n=500)

##ClusterProfiler GO enrichment##
##UPREGULATED GENE LIST###
set.seed(1234)
GO_RvF_up<- enrichGO(names(gene_list_up), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_RvF_up_1 <- clusterProfiler:: simplify(GO_RvF_up, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_RvF_up_1 <- dotplot(GO_RvF_up_1)+ ggtitle(label = "Up_GO_Enrichment: IPR8-R848 v IPR8+flg", subtitle = "All Cells")
dotplt_RvF_up_1 ##Export Plot as PNG 

##########################################################
### Generate list of significantly DOWNREGULATED genes ###
##########################################################

final_RvF_1 <- final_RvF_1[final_RvF_1$p_val_adj < 0.05,]
final_RvF_3 <- final_RvF_1[final_RvF_1$avg_log2FC < 0,]
gene_list_dwn = final_RvF_3$avg_log2FC
names(gene_list_dwn) = final_RvF_3$V1
gene_list_dwn = sort(gene_list_dwn, decreasing = TRUE)
gene_list_dwn = gene_list_dwn[!duplicated(names(gene_list_dwn))] #removes duplicates based on name
head(gene_list_dwn,n=500)

###DOWNREGULATED GENE LIST###
set.seed(1234)
GO_RvF_dwn<- enrichGO(names(gene_list_dwn), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_RvF_dwn_1 <- clusterProfiler:: simplify(GO_RvF_dwn, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_RvF_dwn_1 <- dotplot(GO_RvF_dwn_1)+ ggtitle(label = "Down_GO_Enrichment: IPR8-R848 v IPR8+flg", subtitle = "All Cells")
dotplt_RvF_dwn_1 ##Export Plot as PNG 

##Top 25 most signficantly upregulated genes##
Alltop25_RvF <- final_RvF_genes %>%  
  dplyr::filter(p_val_adj < 0.05) %>% #only significant genes
  dplyr::filter(avg_log2FC > 0) %>% #take upregulated genes
  dplyr::arrange(p_val_adj) %>%  # Sort by p_adj_val 
  slice_head(n = 25) %>%
  ungroup() -> Alltop25_RvF

DotPlot(seurat, features = rownames(Alltop25_RvF), scale = T, dot.scale = 10, group.by = "condition")+
  theme_minimal() + ggtitle(label = "All cells: IPR8-R848 v IPR8+flg", subtitle = "Top 25 most signifcant DEG (Log2FC > 0)") + # Change the theme to minimal and add title 
  theme(plot.background = element_rect(fill = "white")) +  # Set plot background to white
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, margin = margin(t = 10))) +  # Rotate x-axis labels by 45 degrees and adjust horizontal justification
  scale_colour_gradient2(low="black", mid="lightgrey", high="darkgoldenrod")


##Top 25 most significantly downregulated genes##
Allbot25_RvF <- final_RvF_genes %>%  
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr:: filter(avg_log2FC < -0) %>%
  dplyr:: arrange(p_val_adj)%>%
  slice_head(n = 25) %>%
  ungroup() -> Allbot25_RvF


DotPlot(seurat, features = rownames(Allbot25_RvF), scale = T, dot.scale = 10, group.by = "condition")+
  theme_minimal() + ggtitle(label = "All cells: IPR8-R848 v IPR8+flg", subtitle = "Top 25 most signifcant DEG (Log2FC < 0)") +  # Change the theme to minimal
  theme(plot.background = element_rect(fill = "white")) +  # Set plot background to white
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, margin = margin(t = 10))) +  # Rotate x-axis labels by 45 degrees and adjust horizontal justification
  scale_colour_gradient2(low="black", mid="lightgrey", high="darkgoldenrod")