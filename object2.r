#!/usr/bin/env Rscript 
args = commandArgs(trailingOnly=TRUE)

base.path <- args[1]


Read10X <- function( base.path = NULL ){
    if (! dir.exists(base.path )){
      stop("Directory provided does not exist")
    }

    object <- paste0( base.path, "/outs/filtered_feature_bc_matrix/" )
    if (!file.exists( object )){
      stop("input is missing")

require ("suerat")

cellranger.data <- Read10X(base.path)
data <- CreateSeuratObject(raw.data = cellranger.data, min.cells = 3, min.genes = 200, project = "10X_rnaseq")

saveRDS(data, paste0(object, ".rds"))