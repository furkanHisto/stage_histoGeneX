#!/usr/bin/env nextflow

Channel.fromPath(params.bustools).view()
    .set {output_bus}

process bustools_correct{
    input:
    file busfile from output_bus
    file whitelist from params.whitelist

    script:
    """
    bustools correct -w ${whitelist} \
    -o /home/histogenex/bustools/correct/${busfile} \
    ${busfile}
    """
}