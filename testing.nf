#!/usr/bin/env nextflow


// running FASTQC on the samples
/*
fastq_channel = Channel.fromPath('/mnt/storage/Data/10xGenomics_TrainingsRun/*gz')


process fastqc {

    input: 
    
    file fastq from fastq_channel

    """
    fastqc $fastq
    """
}
*/

// executing cellranger count on the files located in 10xGenomics_Trainingsrun


cellranger_channel = Channel.fromPath('/mnt/storage/Data/10xGenomics_TrainingsRun/*.gz')
process cellrangerCount {

    input:
    file cellrangerSample from cellranger_channel

    """
    cellranger count    --id=cellrangerCount --transcriptome=/mnt/storage/Reference/refdata-gex-GRCh38-2020-A --fastqs=/mnt/storage/Data/10xGenomics_TrainingsRun --sample=${cellrangerSample}
    """
}


