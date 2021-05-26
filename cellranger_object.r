#!/usr/bin/env Rscript 

library(Seurat)
library(optparse)

# to gather the input files from the channel. 
#with make option you can write a shebang script on the file. we use this to gather the files emitted from the channel and run the script on it
option_list=list(
  make_option("--sample",default="sample1",action="store",type='character',help="sample"))

opt=parse_args(OptionParser(option_list=option_list))

print(paste0("Sample to process: ", opt$sample))

data <- Read10X(data.dir = ".")
data <-  CreateSeuratObject(counts=data)

saveRDS(data, paste0(opt$sample, "_cellranger_seurat.rds"))