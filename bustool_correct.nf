#!/usr/bin/env nextflow

Channel.fromPath(params.bustools).view()
    .set {output_kbus}


whitelist = Channel.fromPath(params.whitelist)

process bustools_correct{

    publishDir "${params.outdir}/correct/", mode = 'copy'

    input:
    file (busfile) from output_kbus
    file whitelist from whitelist.first()

    script:
    """
    bustools correct -w ${whitelist} \
    -o ./output.bus \
    ${busfile}
    """
}