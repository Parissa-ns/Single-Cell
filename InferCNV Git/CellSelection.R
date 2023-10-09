setwd("D:/scRNASeq/InferCNV")

library(dplyr)
library(Matrix)
library(infercnv)
library(Seurat) 
library(Matrix)
library(ggplot2)
library(dplyr)
library (patchwork) 

 
### CREATING SEURAT object
metadata<-read.csv("D:/scRNASeq/Data/metadata.csv",header=TRUE,row.names=1)
mtx_matrix<-readMM("D:/scRNASeq/Data/matrix.mtx") 
gene_names<-read.table("D:/scRNASeq/Data/genes.tsv", header = FALSE, col.names = 
                           "GeneName", stringsAsFactors = FALSE)$GeneName

rownames(mtx_matrix) <- gene_names
colnames(mtx_matrix) <- rownames(metadata)
# Assign gene names to the matrix
seurat_obj <- CreateSeuratObject(counts= Matrix::Matrix(as.matrix(mtx_matrix),sparse = T), meta.data = metadata)


remove(mtx_matrix)

#InferCNV 
 
celltypesInt<-c("Normal Epithelial", "Cancer Epithelial")  

# Use the celltype_major column to select cells
selected_cells <- which(seurat_obj$celltype_major %in% celltypesInt)

# Subset the Seurat object based on the selected cells
seurat_obj <- seurat_obj[, selected_cells ]

 
############################
###########

# Step 1: Extract the count matrix
count_matrix <- as.data.frame(seurat_obj@assays$RNA@counts ) # Assuming your counts are in the "RNA" assay

# Step 2: Extract the metadata (replace 'metadata_column' with the actual column name you want)
metadata <-  as.data.frame(cbind(colnames(count_matrix), seurat_obj$celltype_major ))
 
# Export the metadata to a CSV file (or any other format you prefer)
write.table(metadata, file = "cellAnnotation.txt",sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Export the count matrix to a CSV file (you can choose a different format if needed)
write.table(count_matrix, file = "count_matrix.txt",sep="\t", row.names = TRUE, col.names=TRUE, quote = FALSE)
 