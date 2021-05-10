#!/usr/bin/env nextflow

Channel     .fromPath(params.outdir3).view()
            .set{bus_sort}

process bus_sort{

    tag "${sort}"

    publishDir "${params.outdir}/kallisto/raw_bus", mode : 'copy'
    
    input:
    file sort from bus_sort


    script:
    """
    bustools sort   -t ${task.cpus} \
                    -o ${sort}.sort.bus \
                    ${sort}
    """
}