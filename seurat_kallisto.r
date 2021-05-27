#!/usr/bin/env Rscript 

#packages
library(BUSpaRse)
library(Seurat)

# to gather the input files from the channel. 
#with make option you can write a shebang script on the file. we use this to gather the files emitted from the channel and run the script on it
args = commandArgs(trailingOnly=TRUE)

sample<- args[1]

#reading file and creating seurat object
res_mat <- read_count_output(sample, name = "tcc", tcc = FALSE)

seu <- CreateSeuratObject(res_mat, min.cells = 3, min.feature=200)

saveRDS(res_mat, paste0(sample, "_kallisto_seurat.rds"))

