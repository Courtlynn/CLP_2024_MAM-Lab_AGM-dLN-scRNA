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

#' ---------- identify DEGs for IPR8+TLRa v PBS ----------- '#
#' ---------- remove Y-linked and polymorphic genes ----------- '#
#' ------ create list of signifcantly up and downregulated genes ----- '#
#' ---------- create GO plot from sig. up/down gene lists ----------- '#


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


################################
## All Cells: IPR8+flg v PBS ##
################################

setwd("~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/DEGs/All Clusters")
Idents(seurat) <-seurat@meta.data$condition
## For the All Clusters analyses below, min.pct is set to 0.01 ##
## Save DEG analysis as csv file after each run ##
FvP <- FindMarkers(seurat, ident.1 = 'Flag', ident.2 = 'PBS', only.pos = FALSE, min.pct = 0.01, logfc.threshold = -Inf, min.cells.feature = 1, min.cells.group = 1)
write.csv(FvP,file='01_AllClusters_FlagvPBS.csv',quote=F)
## Remove Y linked genes ##
filtered_FvP_genes <- FvP[!rownames(FvP) %in% y_hgnc_symbol, ]
##Remove polymorphic HLA-II genes##
final_FvP_genes <- filtered_FvP_genes[!rownames(filtered_FvP_genes) %in% MHCII, ]
write.csv(final_FvP_genes,file='NoY_noMHCII_01_AllClusters_FlagvPBS.csv',quote=F)

######################################################################
## Build the transcript ID to gene ID table from the entrez database##
######################################################################

humanannot<-grch38
humanannot
dfg <- final_FvP_genes
dfg <- as.data.frame(dfg)
dfg$RowNames <- rownames(dfg)
colnames(dfg)[6] <- "RowNames"
df <- cbind(humanannot$entrez, humanannot$symbol)
df <- as.data.frame(df)
dfg
df
final_FvP_1 <- left_join(dfg,df,by=c("RowNames"="V2"))
View(FvP_1)

########################################################
### Generate list of significantly UPERGULATED genes ###
########################################################

final_FvP_1 <- final_FvP_1[final_FvP_1$p_val_adj < 0.05,]
final_FvP_2 <- final_FvP_1[final_FvP_1$avg_log2FC > 0,]
gene_list_up = final_FvP_2$avg_log2FC
names(gene_list_up) = final_FvP_2$V1
gene_list_up = sort(gene_list_up, decreasing = TRUE)
gene_list_up = gene_list_up[!duplicated(names(gene_list_up))] #removes duplicates based on name
head(gene_list_up,n=500)

##ClusterProfiler GO enrichment##
##UPREGULATED GENE LIST###
set.seed(1234)
GO_FvP_up<- enrichGO(names(gene_list_up), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_FvP_up_1 <- clusterProfiler:: simplify(GO_FvP_up, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_FvP_up_1 <- dotplot(GO_FvP_up_1)+ ggtitle(label = "Up_GO_Enrichment: IPR8+flg v PBS", subtitle = "All Cells")
dotplt_FvP_up_1 ##Export Plot as PNG 

##########################################################
### Generate list of significantly DOWNREGULATED genes ###
##########################################################

final_FvP_1 <- final_FvP_1[final_FvP_1$p_val_adj < 0.05,]
final_FvP_3 <- final_FvP_1[final_FvP_1$avg_log2FC < 0,]
gene_list_dwn = final_FvP_3$avg_log2FC
names(gene_list_dwn) = final_FvP_3$V1
gene_list_dwn = sort(gene_list_dwn, decreasing = TRUE)
gene_list_dwn = gene_list_dwn[!duplicated(names(gene_list_dwn))] #removes duplicates based on name
head(gene_list_dwn,n=500)

###DOWNREGULATED GENE LIST###
set.seed(1234)
GO_FvP_dwn<- enrichGO(names(gene_list_dwn), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_FvP_dwn_1 <- clusterProfiler:: simplify(GO_FvP_dwn, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_FvP_dwn_1 <- dotplot(GO_FvP_dwn_1)+ ggtitle(label = "Down_GO_Enrichment: IPR8+flg v PBS", subtitle = "All Cells")
dotplt_FvP_dwn_1 ##Export Plot as PNG 

################################
## All Cells: IPR8-R848 v PBS ##
################################

RvP <- FindMarkers(seurat, ident.1 = 'R848', ident.2 = 'PBS', only.pos = FALSE, min.pct = 0.01, logfc.threshold = -Inf, min.cells.feature = 1, min.cells.group = 1)
write.csv(RvP,file='01_AllClusters_R848vPBS.csv',quote=F)
## Remove Y linked genes ##
filtered_RvP_genes <- RvP[!rownames(RvP) %in% y_hgnc_symbol, ]
##Remove polymorphic HLA-II genes##
final_RvP_genes <- filtered_RvP_genes[!rownames(filtered_RvP_genes) %in% MHCII, ]
write.csv(final_RvP_genes,file='NoY_noMHCII_01_AllClusters_R848vPBS.csv',quote=F)

######################################################################
## Build the transcript ID to gene ID table from the entrez database##
######################################################################

humanannot<-grch38
humanannot
dfg <- final_RvP_genes
dfg <- as.data.frame(dfg)
dfg$RowNames <- rownames(dfg)
colnames(dfg)[6] <- "RowNames"
df <- cbind(humanannot$entrez, humanannot$symbol)
df <- as.data.frame(df)
dfg
df
final_RvP_1 <- left_join(dfg,df,by=c("RowNames"="V2"))
View(RvP_1)

########################################################
### Generate list of significantly UPERGULATED genes ###
########################################################

final_RvP_1 <- final_RvP_1[final_RvP_1$p_val_adj < 0.05,]
final_RvP_2 <- final_RvP_1[final_RvP_1$avg_log2FC > 0,]
gene_list_up = final_RvP_2$avg_log2FC
names(gene_list_up) = final_RvP_2$V1
gene_list_up = sort(gene_list_up, decreasing = TRUE)
gene_list_up = gene_list_up[!duplicated(names(gene_list_up))] #removes duplicates based on name
head(gene_list_up,n=500)

##ClusterProfiler GO enrichment##
##UPREGULATED GENE LIST###
set.seed(1234)
GO_RvP_up<- enrichGO(names(gene_list_up), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_RvP_up_1 <- clusterProfiler:: simplify(GO_RvP_up, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_RvP_up_1 <- dotplot(GO_RvP_up_1)+ ggtitle(label = "Up_GO_Enrichment: IPR8-R848 v PBS", subtitle = "All Cells")
dotplt_RvP_up_1 ##Export Plot as PNG 

##########################################################
### Generate list of significantly DOWNREGULATED genes ###
##########################################################

final_RvP_1 <- final_RvP_1[final_RvP_1$p_val_adj < 0.05,]
final_RvP_3 <- final_RvP_1[final_RvP_1$avg_log2FC < 0,]
gene_list_dwn = final_RvP_3$avg_log2FC
names(gene_list_dwn) = final_RvP_3$V1
gene_list_dwn = sort(gene_list_dwn, decreasing = TRUE)
gene_list_dwn = gene_list_dwn[!duplicated(names(gene_list_dwn))] #removes duplicates based on name
head(gene_list_dwn,n=500)

###DOWNREGULATED GENE LIST###
set.seed(1234)
GO_RvP_dwn<- enrichGO(names(gene_list_dwn), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_RvP_dwn_1 <- clusterProfiler:: simplify(GO_RvP_dwn, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_RvP_dwn_1 <- dotplot(GO_RvP_dwn_1)+ ggtitle(label = "Down_GO_Enrichment: IPR8-R848 v PBS", subtitle = "All Cells")
dotplt_RvP_dwn_1 ##Export Plot as PNG 
