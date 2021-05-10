#!/usr/bin/env nextflow

bus_correct = Channel.fromPath(params.outdir2).view()



whitelist = Channel.fromPath(params.whitelist)

process bustools_correct{
    tag "$bus"

    publishDir "${params.outdir}/kallisto/sort_bus", mode = 'copy'

    input:
    file bus from bus_correct
    file whitelist from whitelist.collect()

    script:
    """
    bustools correct    -w ${whitelist} \
                        -o ${bus}/output.corrected.bus \
                        ${bus}
    """
}