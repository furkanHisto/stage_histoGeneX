#!/usr/bin/env nextflow

cellranger_outs = Channel.fromPath("${params.cellrangerOuts}/*/outs/filtered_feature_bc_matrix/*").view()

process seurat_object {
    tag "${sample}"
    publishDir "${params.outdir}/Seurat",  mode: 'copy' 

    input:
    file sample from cellranger_outs

    output:
    file ".rds" into seurat_paths

    script:
    """
    Rscript /home/histogenex/Pipelines/Nextflow/scRNAseq/seuratObject_cellrangerCount.R --sample ${sample}
    """

}