##Reading in Raw 10X files into RDS##

rm(list=ls())


library(Seurat)
library(Matrix)

############################################################################################
#### set working directory that contains a folders labeled "PBS", "M229", "R848", "Flag" ###
############# in this script, that folder is labeled "Raw Files" ###########################
############################################################################################

##Complete one sample at a time##

##############
#### PBS #####
##############

setwd("~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Raw Files/PBS")
barcodes_path <- "~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Raw Files/PBS/barcodes.tsv.gz"
features_path <- "~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Raw Files/PBS/features.tsv.gz"
matrix_path <- "~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Raw Files/PBS/matrix.mtx.gz"

barcodes <- read.table(barcodes_path, header = FALSE, stringsAsFactors = FALSE)
features <- read.table(features_path, header = FALSE, stringsAsFactors = FALSE)
matrix <- readMM(matrix_path)
#matrix <- as(matrix, "dgCMatrix")


# Ensure feature names are unique
unique_features <- make.unique(features$V2)

# Create the Seurat object
seurat_object <- CreateSeuratObject(counts = matrix, project = "PBS")

# Set row names for the Seurat object
rownames(seurat_object) <- unique_features

# Add cell barcodes as column names for the Seurat object
cell_names <- barcodes$V1
colnames(seurat_object) <- cell_names

# Add the barcodes and condition to the metadata
metadata <- data.frame(barcode = cell_names, condition = "PBS", row.names = cell_names)

# Add the metadata to the Seurat object
seurat_object <- AddMetaData(seurat_object, metadata = metadata)

# Print the Seurat object to confirm
print(seurat_object)

PBS_Seurat<-seurat_object
setwd("~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Objects")
saveRDS(PBS_Seurat,file="PBS_Seurat.RDS")

##############
### M229 #####
##############

setwd("~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Raw Files/M229")
barcodes_path <- "~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Raw Files/M229/barcodes.tsv.gz"
features_path <- "~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Raw Files/M229/features.tsv.gz"
matrix_path <- "~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Raw Files/M229/matrix.mtx.gz"

barcodes <- read.table(barcodes_path, header = FALSE, stringsAsFactors = FALSE)
features <- read.table(features_path, header = FALSE, stringsAsFactors = FALSE)
matrix <- readMM(matrix_path)
#matrix <- as(matrix, "dgCMatrix")


# Ensure feature names are unique
unique_features <- make.unique(features$V2)

# Create the Seurat object
seurat_object <- CreateSeuratObject(counts = matrix, project = "M229")

# Set row names for the Seurat object
rownames(seurat_object) <- unique_features

# Add cell barcodes as column names for the Seurat object
cell_names <- barcodes$V1
colnames(seurat_object) <- cell_names

# Add the barcodes and condition to the metadata
metadata <- data.frame(barcode = cell_names, condition = "M229", row.names = cell_names)

# Add the metadata to the Seurat object
seurat_object <- AddMetaData(seurat_object, metadata = metadata)

# Print the Seurat object to confirm
print(seurat_object)

M229_Seurat<-seurat_object
setwd("~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Objects")
saveRDS(M229_Seurat,file="M229_Seurat.RDS")

##############
### R848 #####
##############

setwd("~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Raw Files/R848")
barcodes_path <- "~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Raw Files/R848/barcodes.tsv.gz"
features_path <- "~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Raw Files/R848/features.tsv.gz"
matrix_path <- "~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Raw Files/R848/matrix.mtx.gz"

barcodes <- read.table(barcodes_path, header = FALSE, stringsAsFactors = FALSE)
features <- read.table(features_path, header = FALSE, stringsAsFactors = FALSE)
matrix <- readMM(matrix_path)
#matrix <- as(matrix, "dgCMatrix")

# Ensure feature names are unique
unique_features <- make.unique(features$V2)

# Create the Seurat object
seurat_object <- CreateSeuratObject(counts = matrix, project = "R848")

# Set row names for the Seurat object
rownames(seurat_object) <- unique_features

# Add cell barcodes as column names for the Seurat object
cell_names <- barcodes$V1
colnames(seurat_object) <- cell_names

# Add the barcodes and condition to the metadata
metadata <- data.frame(barcode = cell_names, condition = "R848", row.names = cell_names)

# Add the metadata to the Seurat object
seurat_object <- AddMetaData(seurat_object, metadata = metadata)

# Print the Seurat object to confirm
print(seurat_object)

R848_Seurat<-seurat_object
setwd("~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Objects")
saveRDS(R848_Seurat,file="R848_Seurat.RDS")

##############
### Flag #####
##############

setwd("~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Raw Files/Flag")
barcodes_path <- "~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Raw Files/Flag/barcodes.tsv.gz"
features_path <- "~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Raw Files/Flag/features.tsv.gz"
matrix_path <- "~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Raw Files/Flag/matrix.mtx.gz"

barcodes <- read.table(barcodes_path, header = FALSE, stringsAsFactors = FALSE)
features <- read.table(features_path, header = FALSE, stringsAsFactors = FALSE)
matrix <- readMM(matrix_path)
#matrix <- as(matrix, "dgCMatrix")


# Ensure feature names are unique
unique_features <- make.unique(features$V2)

# Create the Seurat object
seurat_object <- CreateSeuratObject(counts = matrix, project = "Flag")

# Set row names for the Seurat object
rownames(seurat_object) <- unique_features

# Add cell barcodes as column names for the Seurat object
cell_names <- barcodes$V1
colnames(seurat_object) <- cell_names

# Add the barcodes and condition to the metadata
metadata <- data.frame(barcode = cell_names, condition = "Flag", row.names = cell_names)

# Add the metadata to the Seurat object
seurat_object <- AddMetaData(seurat_object, metadata = metadata)

# Print the Seurat object to confirm
print(seurat_object)

Flag_Seurat<-seurat_object
setwd("~/Desktop/R work/No Dual_scRNA Sequencing Project - CLP and JMG copy/Objects")
saveRDS(Flag_Seurat,file="Flag_Seurat.RDS")
