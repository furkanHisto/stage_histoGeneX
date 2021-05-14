#!/usr/bin/env nextflow

samples_fastqc = Channel.fromPath(params.input)






process fastqc{

    tag "${sample}"

    publishDir "${params.outdir}/fastqc", mode : 'copy' 

    input:
    file sample from samples_fastqc

    output:
    file "*_fastqc*" into fastqc_to_multicq

    script:
    """
    mkdir -p ${params.outdir}/fastqc \

    fastqc  ${sample} \
    
    """  

}

process MultiQC {

    publishDir "${params.outdir}/multiqc", mode : 'copy'

    input:
    file ('fastqc/*') from fastqc_to_multicq.collect().ifEmpty([])

    output:
    file "multiqc_report.html"
    file "multiqc_data"

    """
    multiqc -f . \
    
    """
}

