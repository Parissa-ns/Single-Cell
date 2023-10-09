setwd("D:/scRNASeq/InferCNV")
library(infercnv) 

#create the inferCNV object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="count_matrix.txt",
                                    annotations_file="cellAnnotation.txt",
                                    delim="\t",
                                    gene_order_file="GenomicPositionFile.txt",
                                    ref_group_names=c("Normal Epithelial"),
                                    min_max_counts_per_cell=c(1e3,1e7) )



gc()

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj, 
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics  !!  Mine is 10 x 
                             out_dir="D:/scRNASeq/InferCNV/out",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T, 
                             hclust_method='ward.D2',
                             #k_obs_groups=7, 
                             HMM=T, num_threads = 8
) 
