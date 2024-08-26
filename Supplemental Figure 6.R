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

#' ---------- identify DEGs for IPR8+m229 v PBS ----------- '#
#' ---------- remove Y-linked and polymorphic genes ----------- '#
#' ------ create list of signifcantly up and downregulated genes ----- '#
#' ---------- create GO plot from sig. up/down gene lists ----------- '#
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

##################
### B cells ######
##################
setwd("~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/DEGs/B cells")
seurat$celltype.stim <- paste(seurat$celltype, seurat$condition, sep = "_")
Idents(seurat) <- "celltype.stim"

B_MvP <- FindMarkers(seurat, ident.1 = "B cells_M229", ident.2 = "B cells_PBS", verbose = FALSE)
write.csv(B_MvP,file='B cells_M229vPBS.csv',quote=F)
## Remove Y linked genes ##
filtered_B_MvP_genes <- B_MvP[!rownames(B_MvP) %in% y_hgnc_symbol, ]
##Remove polymorphic HLA-II genes##
final_B_MvP_genes <- filtered_B_MvP_genes[!rownames(filtered_B_MvP_genes) %in% MHCII, ]
write.csv(final_B_MvP_genes,file='NoY_noMHCII_B cells_M229vPBS.csv',quote=F)

######################################################################
## Build the transcript ID to gene ID table from the entrez database##
######################################################################

humanannot<-grch38
humanannot
dfg <- final_B_MvP_genes
dfg <- as.data.frame(dfg)
dfg$RowNames <- rownames(dfg)
colnames(dfg)[6] <- "RowNames"
df <- cbind(humanannot$entrez, humanannot$symbol)
df <- as.data.frame(df)
dfg
df
final_B_MvP_1 <- left_join(dfg,df,by=c("RowNames"="V2"))
View(final_B_MvP_1)

########################################################
### Generate list of significantly UPERGULATED genes ###
########################################################

final_B_MvP_1 <- final_B_MvP_1[final_B_MvP_1$p_val_adj < 0.05,]
final_B_MvP_2 <- final_B_MvP_1[final_B_MvP_1$avg_log2FC > 0,]
gene_list_up = final_B_MvP_2$avg_log2FC
names(gene_list_up) = final_B_MvP_2$V1
gene_list_up = sort(gene_list_up, decreasing = TRUE)
gene_list_up = gene_list_up[!duplicated(names(gene_list_up))] #removes duplicates based on name
head(gene_list_up,n=500)

##ClusterProfiler GO enrichment##
##UPREGULATED GENE LIST###
set.seed(1234)
GO_B_MvP_up<- enrichGO(names(gene_list_up), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_B_MvP_up_1 <- clusterProfiler:: simplify(GO_B_MvP_up, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_B_MvP_up_1 <- dotplot(GO_B_MvP_up_1)+ ggtitle(label = "Up_GO_Enrichment: IPR8+m229 v PBS", subtitle = "B cells")
dotplt_B_MvP_up_1 ##Export Plot as PNG 

##########################################################
### Generate list of significantly DOWNREGULATED genes ###
##########################################################

final_B_MvP_1 <- final_B_MvP_1[final_B_MvP_1$p_val_adj < 0.05,]
final_B_MvP_3 <- final_B_MvP_1[final_B_MvP_1$avg_log2FC < 0,]
gene_list_dwn = final_B_MvP_3$avg_log2FC
names(gene_list_dwn) = final_B_MvP_3$V1
gene_list_dwn = sort(gene_list_dwn, decreasing = TRUE)
gene_list_dwn = gene_list_dwn[!duplicated(names(gene_list_dwn))] #removes duplicates based on name
head(gene_list_dwn,n=500)

###DOWNREGULATED GENE LIST###
set.seed(1234)
GO_B_MvP_dwn<- enrichGO(names(gene_list_dwn), 'org.Hs.eg.db', ont="BP", keyType = "ENTREZID", pvalueCutoff=0.01) 
GO_B_MvP_dwn_1 <- clusterProfiler:: simplify(GO_ISG_MvP_dwn, cutoff=0.7, by="p.adjust", select_fun=min)
dotplt_B_MvP_dwn_1 <- dotplot(GO_B_MvP_dwn_1)+ ggtitle(label = "Down_GO_Enrichment: IPR8+m229 v PBS", subtitle = "B cells")
dotplt_B_MvP_dwn_1 ##Export Plot as PNG

