#!/usr/bin/env Rscript 

library(Seurat)

# to gather the input files from the channel. 
#with make option you can write a shebang script on the file. we use this to gather the files emitted from the channel and run the script on it
args = commandArgs(trailingOnly=TRUE)

print(args)

sample<- args[1]

#loading  files

data <- Read10X(data.dir = ".")
data <-  CreateSeuratObject(counts=data, min.cells = 3, min.features = 200, project = paste0(sample, "_cellranger"))



#save seurat obj into a .rds file
saveRDS(data, paste0(sample, "_cellranger_seurat.rds"))