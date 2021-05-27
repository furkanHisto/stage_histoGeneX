#!/usr/bin/env nextflow


//creating channel to get the kallisto index file, the whitelist and the GTF file of the v32 reference genome
kallisto_index = Channel.fromPath(params.kallisto_index)
                        .ifEmpty { exit 1, "Index not found: ${params.kallisto_index}" }

whitelist = Channel.fromPath(params.whitelist)

gtf_gene_map = Channel.fromPath(params.gtf_gene_map)

//samples_fastqc = Channel.fromPath("${params.input_kallisto}")




//channel to read the sample files
Channel .fromFilePairs(params.input_cellranger2).view()                                                     // the files emited are collected in pairs into a tuple. 
        .ifEmpty {exit 1, "params.input_paths was empty - no input files supplied" }
        .map {tag, pair ->                                                              // each sample occurs twice. this needs to change to 1.  a tag opperator is used to group the samples
                subtags = (tag =~ /sample\d{1}/);                                       // first you write a pattern to match folowed by \d to specify a wild card.
                tuple ( subtags[0] ,pair )                                              //here you specify how the tuple is ordered. you put the matching string tag as the first value in the tuple
        }
        .groupTuple(by:[0])                                                        // you then sort by tuple. you can specify by which tuple value you want to sort. because we specified in the earlier step that sample\d is the first one, it will sort by sample by default.
        .map{sample, file ->
                tuple(sample, file[0],file[1])
        }
        .set {read_files_kallisto}
        



// process for the kallisto bus command. this will create a matrix.ec file, a transcript.txt file and a output.bus file. it counts all the transcripts.

process kallisto_bus{

        label 'mid_memory'
        
        tag "${name}"

        publishDir "${params.outdir}/kallisto/raw_bus2", mode: 'copy'

        input:
        tuple val (name), file (read1), file (read2) from read_files_kallisto
        file index from kallisto_index.collect()
        


        output:
        file "${name}_bus_output" into kallisto_bus_to_correct
        file "${name}_kallisto.log" into kallisto_logs

        

         """
         kallisto bus  -i $index \
                       --output=${name}_bus_output/ \
                       -x '10xv2' -t ${task.cpus} \
                       ${read1} ${read2} 2>&1 | tee ${name}_kallisto.log
        """
}

// process for bustools correct.  this will filter the counts to a whitelist so that only the unique counts are counted.


process bustools_correct{

        label 'mid_memory'

        tag "$bus"

        publishDir "${params.outdir}/kallisto/raw_bus2", mode : 'copy'

        input:
        file bus from kallisto_bus_to_correct
        file whitelist from whitelist.collect()

        output:
        file bus into corrected_to_sort

        script:
        """
        bustools correct    -w ${whitelist} \
                                -o ${bus}/output.corrected.bus \
                                ${bus}/output.bus
        """
}


// proces for bustools sort command. this will sort the bus files so that the downwards processes can be executed faster.
process bustools_sort{

        label 'mid_memory'

        tag "${sort}"

        publishDir "${params.outdir}/kallisto/raw_bus2", mode : 'copy'
        
        input:
        file sort from corrected_to_sort

        output:
        file sort into sort_to_count
        file sort into sort_to_inspect


        script:
        """
        bustools sort       -t ${task.cpus} \
                                -o ${sort}/output.corrected.sorted.bus \
                                -m ${task.memory.toGiga()} \
                                ${sort}/output.corrected.bus
        """
}

// process of the Bustools inspect command. this will create a sumary file of the bus files in json format.
process bustools_inspect{

        tag "${inspect}"

        publishDir "${params.outdir}/kallisto/inspect2", mode : 'copy'
        
        input:
        file inspect from sort_to_inspect

        output:
        file "${inspect}.json"

        script:
        """
        bustools inspect    -o ${inspect}.json \
                                ${inspect}/output.corrected.sorted.bus
        """
}




// process to make a gene map from the GTF file from the reference genome v32.
process gene_map{

        tag "$gtf"

        publishDir "${params.outdir}/kallisto/gene_map", mode : 'copy'

        input:
        file gtf from gtf_gene_map

        output:
        file "transcript_to_gene.txt" into kallisto_gene_map
        
        shell:
        '''
        zless -S !{gtf} | grep -v "#" | \
        awk '$3 =="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$4"\t"$2}' | sort \
        | uniq |  sed 's/\"//g' | tee transcript_to_gene.txt
        '''
}

process bustools_count {

        tag "$bus"

        publishDir "${params.outdir}/kallisto/count2", mode : 'copy'

        input:
        file bus from sort_to_count
        file gene_map from kallisto_gene_map.collect()

        output:
        file("${bus}_count")  into kallisto_for_seurat

        script:
        """
        mkdir -p ${bus}_count
        mkdir -p ${bus}_genecount
        bustools count  -o ${bus}_count/tcc \
                        -g ${gene_map} \
                        -e ${bus}/matrix.ec \
                        -t ${bus}/transcripts.txt \
                        ${bus}/output.corrected.sorted.bus

        bustools count  -o ${bus}_count/tcc \
                        -g ${gene_map} \
                        -e ${bus}/matrix.ec \
                        -t ${bus}/transcripts.txt \
                        --genecounts \
                        ${bus}/output.corrected.sorted.bus
        """
        
}


process kallisto_seurat_object {
    tag "${bus}"
    publishDir "${params.outdir}/Seurat2",  mode: 'copy' 

    input:
    file (bus) from kallisto_for_seurat

    output:
    file "*.rds"

    script:
    """
    Rscript /home/histogenex/Pipelines/Nextflow/scRNAseq/seurat_kallisto.r   ${bus}
    """

}
