#!/usr/bin/env Rscript 

#install.packages("Seurat", repos = "https://satijalab.org/seurat")
library(Seurat)
#if (!require(optparse)) install.packages("optparse",repos="http://cran.rstudio.com/",quiet=TRUE)
#library(optparse)

####
# to gather the input files from the channel. 
#with make option you can write a shebang script on the file. we use this to gather the files emitted from the channel and run the script on it
#option_list=list(
#  make_option("--sample",default="sample1",action="store",type='character',help="sample"))

#opt=parse_args(OptionParser(option_list=option_list))
#rm(option_list)


####
#trying to get data second option
library (argparser)

ParseArguments <- function() {
  p <- arg_parser('Run Seurat')
  p <- add_argument(p, 'feature-matrix',
                      help='folder or h5 containing feature matrix')
  return(parse_args(p))
  p <- add_argument(p, '--aggregation',
                      help='csv containing metadata, output by Cell Ranger')
  p <- add_argument(p, '--output-dir', default='.',
                      help='output directory for plots and tables')
}

arg <- ParseArguments()

# create output directory
dir.create(arg$output_dir)

#loading data
data_dir <- Read10X(data.dir = arg$feature_matrix)
#dim(data_dir)
#create seurat object of the data
seurat <-  CreateSeuratObject(counts = data_dir)

#testing multpile seurat commands on object
#data <- NormalizeData(object = data)
#data <- FindVariableFeatures(object = data)
#data <- ScaleData(object = data)
#data <- RunPCA(object = data)
#data <- FindNeighbors(object = data)
#data <- FindClusters(object = data)
#data <- RunTSNE(object = data)


# Get cell and feature names, and total numbers
#colnames(x = data)
#rownames(x = data)
#ncol(x = data)
#nrow(x = data)
  
# View metadata data frame, stored in object@meta.data
#data[[]]

#saving object into a file
saveRDS(seurat, file = file.path(arg$output_dir, 'seurat.rds'))