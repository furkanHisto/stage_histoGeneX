#!/usr/bin/env nextflow

// gtf_gene_map = Channel.fromPath(params.gtf_gene_map)
sort_to_count = Channel.fromPath(params.outdir2).view()
kallisto_gene_map = Channel.fromPath(params.gene_map)


process to make a gene map from the GTF file from the reference genome v32
process gene_map{

    tag "$gtf"

    publishDir "${params.outdir}/kallisto/gene_map", mode : 'copy'

    input:
    file gtf from gtf_gene_map

    output:
    file "transcript_to_gene.txt" into kallisto_gene_map
    
    shell:
    '''
    zless -S /mnt/storage/Reference/GENCODE/v32/GTF/gencode.v32.primary_assembly.annotation.gtf | grep -v "#" | \
    awk '$3 =="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$4"\t"$2}' | sort \
    | uniq |  sed 's/\"//g' | tee transcript_to_gene.txt
    '''
}


}