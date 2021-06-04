#!/usr/bin/env Rscript 

#packages
library(BUSpaRse)
library(Seurat)

# to gather the input files from the channel. 
#with make option you can write a shebang script on the file. we use this to gather the files emitted from the channel and run the script on it
args = commandArgs(trailingOnly=TRUE)

sample<- args[1]

#reading file 
res_mat <- read_count_output(".", name = "gene", tcc = FALSE)
#rename "sample" to create better filenames
renamed <- gsub( "_bus_output", "", sample)
#creating seurat object
seu <- CreateSeuratObject(res_mat, min.cells = 3, min.feature=200, project = renamed)


#save seurat obj into a .rds file
saveRDS(seu, paste0(renamed, "_seurat.rds"))

