#!/usr/bin/env nextflow

cellranger_outs = Channel.fromPath("${params.cellrangerOuts}/*", type: 'dir')

process seurat_object {
    tag "${sample}"
    publishDir "${params.outdir}/Seurat",  mode: 'copy' 

    input:
    val sample, file cellranger_barcodes, cellranger_features, cellranger_matrix from cellranger_outs

    output:
    file "*.rds"

    script:
    """
    Rscript /home/histogenex/Pipelines/Nextflow/scRNAseq/cellranger_object.r --sample ${sample}
    """

}