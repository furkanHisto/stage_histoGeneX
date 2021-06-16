# stage_histoGeneX

A pipeline created using nextflow to preprocess 10x Genomic single-cell RNA-seq data into Seurat object files.

2 workflows are done to preprocess the data, cell ranger and kallisto. The output of these workflows are converted into a Seurat object using a custom-made R script.
