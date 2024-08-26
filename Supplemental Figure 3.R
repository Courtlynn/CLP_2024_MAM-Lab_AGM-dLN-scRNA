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

#' ----- DEGs for IPR8+m229 and IPR8+TLRa v PBS ---- '#
#' ---------- remove Y-linked and polymorphic genes ----------- '#
#' ------ create list of signifcantly up and downregulated genes ----- '#
#' ------- create Dotplot of GO Enrichment Results -------- '#

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



##################################
### APCs: IPR8+m229 v PBS ########
##################################
setwd("~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/DEGs/APCs")
Idents(seurat) <- "celltype.stim"
APC_MvP <- FindMarkers(seurat, ident.1 = c("Dendritic cells_M229", "Monocytes/Macs_M229"), ident.2 = c("Dendritic cells_PBS", "Monocytes/Macs_PBS"), verbose = FALSE)
write.csv(APC_MvP,file='APCs_M229vPBS.csv',quote=F)
## Remove Y linked genes ##
filtered_APC_MvP_genes <- APC_MvP[!rownames(APC_MvP) %in% y_hgnc_symbol, ]
##Remove polymorphic HLA-II genes##
final_APC_MvP_genes <- filtered_APC_MvP_genes[!rownames(filtered_APC_MvP_genes) %in% MHCII, ]
write.csv(final_APC_MvP_genes,file='NoY_noMHCII_APC_M229vPBS.csv',quote=F)

######################################################################
## Build the transcript ID to gene ID table from the entrez database##
######################################################################

humanannot<-grch38
humanannot
dfg <- final_APC_MvP_genes
dfg <- as.data.frame(dfg)
dfg$RowNames <- rownames(dfg)
colnames(dfg)[6] <- "RowNames"
df <- cbind(humanannot$entrez, humanannot$symbol)
df <- as.data.frame(df)
dfg
df
final_APC_MvP_1 <- left_join(dfg,df,by=c("RowNames"="V2"))
View(MvP_1)

########################################################
### Generate list of significantly UPERGULATED genes ###
########################################################

final_APC_MvP_1 <- final_APC_MvP_1[final_APC_MvP_1$p_val_adj < 0.05,]
final_APC_MvP_2 <- final_APC_MvP_1[final_APC_MvP_1$avg_log2FC > 0,]
gene_list_up = final_APC_MvP_2$avg_log2FC
names(gene_list_up) = final_APC_MvP_2$V1
gene_list_up = sort(gene_list_up, decreasing = TRUE)
gene_list_up = gene_list_up[!duplicated(names(gene_list_up))] #removes duplicates based on name
#gene_list = gene_list[!duplicated(gene_list)] #removes duplicates based on value - should we include this???
head(gene_list_up,n=500)

##ClusterProfiler GO enrichment##
##UPREGULATED GENE LIST###
set.seed(1234)
GO_APC_MvP_up<- enrichGO(names(gene_list_up), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.05) 
GO_APC_MvP_up_1 <- clusterProfiler:: simplify(GO_APC_MvP_up, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_APC_MvP_up_1 <- dotplot(GO_APC_MvP_up_1)+ ggtitle(label = "Up_GO_Enrichment: IPR8+m229 v PBS", subtitle = "APCs")
dotplt_APC_MvP_up_1 ##Export plot as PNG

##########################################################
### Generate list of significantly DOWNREGULATED genes ###
##########################################################

final_APC_MvP_1 <- final_APC_MvP_1[final_APC_MvP_1$p_val_adj < 0.05,]
final_APC_MvP_3 <- final_APC_MvP_1[final_APC_MvP_1$avg_log2FC < 0,]
gene_list_dwn = final_APC_MvP_3$avg_log2FC
names(gene_list_dwn) = final_APC_MvP_3$V1
gene_list_dwn = sort(gene_list_dwn, decreasing = TRUE)
gene_list_dwn = gene_list_dwn[!duplicated(names(gene_list_dwn))] #removes duplicates based on name
#gene_list = gene_list[!duplicated(gene_list)] #removes duplicates based on value - should we include this???
head(gene_list_dwn,n=500)

###DOWNREGULATED GENE LIST###
set.seed(1234)
GO_APC_MvP_dwn<- enrichGO(names(gene_list_dwn), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.05) 
GO_APC_MvP_dwn_1 <- clusterProfiler:: simplify(GO_APC_MvP_dwn, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_APC_MvP_dwn_1 <- dotplot(GO_APC_MvP_dwn_1)+ ggtitle(label = "Down_GO_Enrichment: IPR8+m229 v PBS", subtitle = "APCs")
dotplt_APC_MvP_dwn_1 ##Export plot as PNG

##################################
### APCs: IPR8+flg v PBS #########
##################################
APC_FvP <- FindMarkers(seurat, ident.1 = c("Dendritic cells_Flag", "Monocytes/Macs_Flag"), ident.2 = c("Dendritic cells_PBS", "Monocytes/Macs_PBS"), verbose = FALSE)
write.csv(APC_FvP,file='APCs_FlagvPBS.csv',quote=F)
## Remove Y linked genes ##
filtered_APC_FvP_genes <- APC_FvP[!rownames(APC_FvP) %in% y_hgnc_symbol, ]
##Remove polymorphic HLA-II genes##
final_APC_FvP_genes <- filtered_APC_FvP_genes[!rownames(filtered_APC_FvP_genes) %in% MHCII, ]
write.csv(final_APC_FvP_genes,file='NoY_noMHCII_APC_FlagvPBS.csv',quote=F)

######################################################################
## Build the transcript ID to gene ID table from the entrez database##
######################################################################

humanannot<-grch38
humanannot
dfg <- final_APC_FvP_genes
dfg <- as.data.frame(dfg)
dfg$RowNames <- rownames(dfg)
colnames(dfg)[6] <- "RowNames"
df <- cbind(humanannot$entrez, humanannot$symbol)
df <- as.data.frame(df)
dfg
df
final_APC_FvP_1 <- left_join(dfg,df,by=c("RowNames"="V2"))
View(final_APC_FvP_1)

########################################################
### Generate list of significantly UPERGULATED genes ###
########################################################

final_APC_FvP_1 <- final_APC_FvP_1[final_APC_FvP_1$p_val_adj < 0.05,]
final_APC_FvP_2 <- final_APC_FvP_1[final_APC_FvP_1$avg_log2FC > 0,]
gene_list_up = final_APC_FvP_2$avg_log2FC
names(gene_list_up) = final_APC_FvP_2$V1
gene_list_up = sort(gene_list_up, decreasing = TRUE)
gene_list_up = gene_list_up[!duplicated(names(gene_list_up))] #removes duplicates based on name
#gene_list = gene_list[!duplicated(gene_list)] #removes duplicates based on value - should we include this???
head(gene_list_up,n=500)

##ClusterProfiler GO enrichment##
##UPREGULATED GENE LIST###
set.seed(1234)
GO_APC_FvP_up<- enrichGO(names(gene_list_up), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.05) 
GO_APC_FvP_up_1 <- clusterProfiler:: simplify(GO_APC_FvP_up, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_APC_FvP_up_1 <- dotplot(GO_APC_FvP_up_1)+ ggtitle(label = "Up_GO_Enrichment: IPR8+flg v PBS", subtitle = "APCs")
dotplt_APC_FvP_up_1 ##Export plot as PNG

##########################################################
### Generate list of significantly DOWNREGULATED genes ###
##########################################################

final_APC_FvP_1 <- final_APC_FvP_1[final_APC_FvP_1$p_val_adj < 0.05,]
final_APC_FvP_3 <- final_APC_FvP_1[final_APC_FvP_1$avg_log2FC < 0,]
gene_list_dwn = final_APC_FvP_3$avg_log2FC
names(gene_list_dwn) = final_APC_FvP_3$V1
gene_list_dwn = sort(gene_list_dwn, decreasing = TRUE)
gene_list_dwn = gene_list_dwn[!duplicated(names(gene_list_dwn))] #removes duplicates based on name
#gene_list = gene_list[!duplicated(gene_list)] #removes duplicates based on value - should we include this???
head(gene_list_dwn,n=500)

###DOWNREGULATED GENE LIST###
set.seed(1234)
GO_APC_FvP_dwn<- enrichGO(names(gene_list_dwn), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.05) 
GO_APC_FvP_dwn_1 <- clusterProfiler:: simplify(GO_APC_FvP_dwn, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_APC_FvP_dwn_1 <- dotplot(GO_APC_FvP_dwn_1)+ ggtitle(label = "Down_GO_Enrichment: IPR8+flg v PBS", subtitle = "APCs")
dotplt_APC_FvP_dwn_1 ##Export plot as PNG


##################################
### APCs: IPR8-R848 v PBS ########
##################################
APC_RvP <- FindMarkers(seurat, ident.1 = c("Dendritic cells_R848", "Monocytes/Macs_R848"), ident.2 = c("Dendritic cells_PBS", "Monocytes/Macs_PBS"), verbose = FALSE)
write.csv(APC_RvP,file='APCs_R848vPBS.csv',quote=F)
## Remove Y linked genes ##
filtered_APC_RvP_genes <- APC_RvP[!rownames(APC_RvP) %in% y_hgnc_symbol, ]
##Remove polymorphic HLA-II genes##
final_APC_RvP_genes <- filtered_APC_RvP_genes[!rownames(filtered_APC_RvP_genes) %in% MHCII, ]
write.csv(final_APC_RvP_genes,file='NoY_noMHCII_APC_R848vPBS.csv',quote=F)

######################################################################
## Build the transcript ID to gene ID table from the entrez database##
######################################################################

humanannot<-grch38
humanannot
dfg <- final_APC_RvP_genes
dfg <- as.data.frame(dfg)
dfg$RowNames <- rownames(dfg)
colnames(dfg)[6] <- "RowNames"
df <- cbind(humanannot$entrez, humanannot$symbol)
df <- as.data.frame(df)
dfg
df
final_APC_RvP_1 <- left_join(dfg,df,by=c("RowNames"="V2"))
View(final_APC_RvP_1)

########################################################
### Generate list of significantly UPERGULATED genes ###
########################################################

final_APC_RvP_1 <- final_APC_RvP_1[final_APC_RvP_1$p_val_adj < 0.05,]
final_APC_RvP_2 <- final_APC_RvP_1[final_APC_RvP_1$avg_log2FC > 0,]
gene_list_up = final_APC_RvP_2$avg_log2FC
names(gene_list_up) = final_APC_RvP_2$V1
gene_list_up = sort(gene_list_up, decreasing = TRUE)
gene_list_up = gene_list_up[!duplicated(names(gene_list_up))] #removes duplicates based on name
#gene_list = gene_list[!duplicated(gene_list)] #removes duplicates based on value - should we include this???
head(gene_list_up,n=500)

##ClusterProfiler GO enrichment##
##UPREGULATED GENE LIST###
set.seed(1234)
GO_APC_RvP_up<- enrichGO(names(gene_list_up), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.05) 
GO_APC_RvP_up_1 <- clusterProfiler:: simplify(GO_APC_RvP_up, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_APC_RvP_up_1 <- dotplot(GO_APC_RvP_up_1)+ ggtitle(label = "Up_GO_Enrichment: IPR8-R848 v PBS", subtitle = "APCs")
dotplt_APC_RvP_up_1 ##Export plot as PNG

##########################################################
### Generate list of significantly DOWNREGULATED genes ###
##########################################################

final_APC_RvP_1 <- final_APC_RvP_1[final_APC_RvP_1$p_val_adj < 0.05,]
final_APC_RvP_3 <- final_APC_RvP_1[final_APC_RvP_1$avg_log2FC < 0,]
gene_list_dwn = final_APC_RvP_3$avg_log2FC
names(gene_list_dwn) = final_APC_RvP_3$V1
gene_list_dwn = sort(gene_list_dwn, decreasing = TRUE)
gene_list_dwn = gene_list_dwn[!duplicated(names(gene_list_dwn))] #removes duplicates based on name
#gene_list = gene_list[!duplicated(gene_list)] #removes duplicates based on value - should we include this???
head(gene_list_dwn,n=500)

###DOWNREGULATED GENE LIST###
set.seed(1234)
GO_APC_RvP_dwn<- enrichGO(names(gene_list_dwn), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.05) 
GO_APC_RvP_dwn_1 <- clusterProfiler:: simplify(GO_APC_RvP_dwn, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_APC_RvP_dwn_1 <- dotplot(GO_APC_RvP_dwn_1)+ ggtitle(label = "Down_GO_Enrichment: IPR8-R848 v PBS", subtitle = "APCs")
dotplt_APC_RvP_dwn_1 ##Export plot as PNG

