#!/usr/bin/env Rscript 

#packages

library(Seurat)

# to gather the input files from the channel. 
#with make option you can write a shebang script on the file. we use this to gather the files emitted from the channel and run the script on it
args = commandArgs(trailingOnly=TRUE)
sample<- args[1]  


#adding second argument to specify which process it is
process_specified <- args[2]

#create the right object specified by the process


write.table(process_specified, "current_process.txt")

if (process_specified == "cellranger") {
  data <- Read10X(data.dir = ".")
  data <-  CreateSeuratObject(counts=data, min.cells = 3, min.features = 200, project = paste0(sample, "_cellranger"))
  saveRDS(data, paste0(sample, "_cellranger_seurat.rds"))
} else if (process_specified == "kallisto_tcc")  {
  require(BUSpaRse)
  res_mat <- read_count_output(".", name = "tcc", tcc = TRUE)
  renamed <- gsub( "_bus_output", "", sample)
  seu <- CreateSeuratObject(res_mat, min.cells = 3, min.feature=200, project = renamed)
  saveRDS(seu, paste0(renamed, "_kallisto_seurat.rds"))
} else if (process_specified == "kallisto_gene")  {
  require(BUSpaRse)
  res_mat <- read_count_output(".", name = "gene", tcc = FALSE)
  renamed <- gsub( "_bus_output", "", sample) 
  seu <- CreateSeuratObject(res_mat, min.cells = 3, min.feature=200, project = renamed)
  saveRDS(seu, paste0(renamed, "_kallisto_seurat.rds"))
}



#option2
# datas <- function(fileDir =NULL){
#     if (! file.exists("*.gz")){
#         data <- Read10X(data.dir = ".")
#         data <-  CreateSeuratObject(counts=data, min.cells = 3, min.features = 200, project = paste0(sample, "_cellranger") 
#         return(data)
#     } else if (! file.exists("tcc*")){
#         #reading file 
#         res_mat <- read_count_output(".", name = "tcc", tcc = TRUE)
#         #rename "sample" to create better filenames
#         renamed <- gsub( "_bus_output", "", sample)
#         #creating seurat object
#         seu <- CreateSeuratObject(res_mat, min.cells = 3, min.feature=200, project = renamed)
#         return (seu)
#     }
#     else (! file.exists('gene*')){
#         #reading file 
#         res_mat <- read_count_output(".", name = "gene", tcc = FALSE)
#         #rename "sample" to create better filenames
#         renamed <- gsub( "_bus_output", "", sample)
#         #creating seurat object
#         seu <- CreateSeuratObject(res_mat, min.cells = 3, min.feature=200, project = renamed)
#     }
# }





#if statements to create seurat object dependant on which file its read

# for (file in fileDir)
#   if (grep("*.gz", file)){
#     #reading files
#     seu1 <- Read10X(data.dir = ".")
#     seu1 <-  CreateSeuratObject(counts=seu, min.cells = 3, min.features = 200, project = paste0(sample, "_cellranger"))
#     #save seurat obj into a .rds file
#     saveRDS(seu1, paste0(sample, "_cellranger_seurat.rds"))                        
#   } else if (grep("tcc*", file)){
#     #reading file 
#     res_mat <- read_count_output(".", name = "tcc", tcc = TRUE)
#     #rename "sample" to create better filenames
#     renamed <- gsub( "_bus_output", "", sample)
#     #creating seurat object
#     seu2 <- CreateSeuratObject(res_mat, min.cells = 3, min.feature=200, project = renamed)
#     #save seurat obj into a .rds file
#     saveRDS(seu2, paste0(renamed, "_kallisto_seurat.rds"))
#   } else (grep("gene*", file)){
#     #reading file 
#     res_mat2 <- read_count_output(".", name = "gene", tcc = FALSE)
#     #rename "sample" to create better filenames
#     renamed2 <- gsub( "_bus_output", "", sample)
#     #creating seurat object
#     seu3 <- CreateSeuratObject(res_mat, min.cells = 3, min.feature=200, project = renamed2)
#     #save seurat obj into a .rds file
#     saveRDS(seu3, paste0(renamed2, "_seurat.rds"))
#   }
