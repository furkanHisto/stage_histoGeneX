#!/usr/bin/env nextflow

// samples_fastqc = Channel.fromPath(params.input)

sample_multiqc = Channel.fromPath(params.sample_multiqc).view()


process MultiQC {

    publishDir "${params.outdir}/multiqc", mode : 'copy'

    input:
    file sample from sample_multiqc    

    output:
    file "${sample}"

    """
    multiqc . \
    
    """
}

// process fastqc{

//     tag "${sample}"

//     publishDir "${params.outdir}/fastqc", mode : 'copy' 

//     input:
//     file sample from samples_fastqc

//     output:
//     file "${sample}"

//     script:
//     """
//     mkdir -p ${params.outdir}/fastqc \

//     fastqc --extract \
//             -o ${params.outdir}/fastqc \
//             ${sample} \
    
//     """  

// }

