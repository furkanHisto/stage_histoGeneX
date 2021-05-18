install.packages("Seurat")
library(Seurat)
install.packages("optparse")
library(optparse)

# to gather the input files from the channel. 
#with make option you can write a shebang script on the file. we use this to gather the files emitted from the channel and run the script on it
option_list=list(
  make_option("--sample",default="sample1",action="store",type='character',help="sample"))

opt=parse_args(OptionParser(option_list=option_list))
rm(option_list)

#loading data
data_dir <- Read10X(data.dir = "opt/")
dim(data_dir)
#create seurat object of the data
data <-  CreateSeuratObject(counts = data_dir)

#testing multpile seurat commands on object
data <- NormalizeData(object = data)
data <- FindVariableFeatures(object = data)
data <- ScaleData(object = data)
data <- RunPCA(object = data)
data <- FindNeighbors(object = data)
data <- FindClusters(object = data)
data <- RunTSNE(object = data)


# Get cell and feature names, and total numbers
colnames(x = data)
rownames(x = data)
ncol(x = data)
nrow(x = data)
  
# View metadata data frame, stored in object@meta.data
data[[]]

#saving object into a file
saveRDS(data, file = "sample1.rds")
