#!/usr/bin/env nextflow

// channel to get the transcriptome file for use in cellranger count
transcriptome_reference = Channel.fromPath(params.transcriptome)
                                 .ifEmpty { exit 1, "Reference not found: ${params.transcriptome}" }


// cahnnel to get the path of the samples. cellranger count needs the path of the samples to work
fastqdir = Channel.fromPath(params.input_cellranger)
                  .ifEmpty { exit 1, "Samples not found: ${params.input_cellranger}" }


//channel to get the sample files for cellranger count
sampleChannel = Channel.fromPath("${params.input_cellranger}/*.gz")
    .ifEmpty { exit 1, "Samples not found: ${params.input_cellranger}" }
    .map { Path path -> 
        path.toFile()
            .getSimpleName()                                       // this returns the name without the extension
            .replaceAll(/_S[1-4]_L00[1-2]_[IR][1-2]_001/, '')      // this replaces the first string with something else. we need to change a part of the file name so the cellranger count can recognize the file name
    }
    .unique()  //this makes it so that there are no duplicate files



//channel to get the kallisto index file, the whitelist and the GTF file of the v32 reference genome
kallisto_index = Channel.fromPath(params.kallisto_index)
                        .ifEmpty { exit 1, "Index not found: ${params.kallisto_index}" }

whitelist = Channel.fromPath(params.whitelist)

gtf_gene_map = Channel.fromPath(params.gtf_gene_map)

samples_fastqc = Channel.fromPath("${params.input_kallisto}")




//channel to read the sample files for kallisto bus. 
Channel .fromFilePairs(params.input_kallisto)                                                     // the files emited are collected in pairs into a tuple. 
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

        publishDir "${params.outdir}/kallisto/raw_bus", mode: 'copy'

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

        publishDir "${params.outdir}/kallisto/raw_bus", mode : 'copy'

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

        publishDir "${params.outdir}/kallisto/raw_bus", mode : 'copy'
        
        input:
        file sort from corrected_to_sort

        output:
        file sort into sort_to_count
        file sort into sort_to_inspect


        script:
        """
        bustools sort   -t ${task.cpus} \
                        -o ${sort}/output.corrected.sorted.bus \
                        -m ${task.memory.toGiga()} \
                        ${sort}/output.corrected.bus
        """
}

// process of the Bustools inspect command. this will create a sumary file of the bus files in json format.
process bustools_inspect{

        tag "${inspect}"

        publishDir "${params.outdir}/kallisto/inspect", mode : 'copy'
        
        input:
        file inspect from sort_to_inspect

        output:
        file "${inspect}.json"

        script:
        """
        bustools inspect        -o ${inspect}.json \
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

        publishDir "${params.outdir}/kallisto/count", mode : 'copy'

        input:
        file bus from sort_to_count
        file gene_map from kallisto_gene_map.collect()

        output:
        file "${bus}_count"

        script:
        """
        mkdir -p ${bus}_count
        bustools count  -o ${bus}_count/tcc \
                        -g ${gene_map} \
                        -e ${bus}/matrix.ec \
                        -t ${bus}/transcripts.txt \
                        ${bus}/output.corrected.sorted.bus
        """
        
}


// process for creating a FASTQC and a MultiQC report of the samples.
process fastqc{

        tag "${sample}"

        publishDir "${params.outdir}/fastqc", mode : 'copy' 

        input:
        file sample from samples_fastqc

        output:
        file "*_fastqc*" into fastqc_to_multicq

        script:
        """
        fastqc  ${sample} \
        
        """  

}


process MultiQC {

        publishDir "${params.outdir}/multiqc", mode : 'copy'

        input:
        file ('fastqc/*') from fastqc_to_multicq.collect().ifEmpty([])
        // file ('kallisto.log') from kallisto_logs.collect().ifEmpty([])

        output:
        file "multiqc_report.html"
        file "multiqc_data"

        """
        multiqc -f  . 
        
        """
}


// process for cellranger count.

process CellRangerCount {

  publishDir "${params.outdir}/cellrangercount",  mode: 'copy'       // this is added to specify the directory of the output of the run.

  label "mid_memory"     // this looks in the hgx.config file for this label. this label sets the configuration of all the processes to set values.
  tag "$sample"         // this is added to check wich sample is currently running when executing nextflow.

  input:
  val sample from sampleChannel      // we take the values generated from the samplechannel. these are the file names 'GEM_sample*' 
  file fastq_dir from fastqdir.first()              // we add .first because that channel only outputs 1 file while the samplechannel outputs 4. when there is no .first() then the run will stop at 1 file because the fastq_dir channel only outputs 1 file. 
  file tx_reference from transcriptome_reference.first()    // .first takes the first output file and reuses it. so the only output of these channels is reuses 4 times to match the amount of samplechannel output.
  
  output:
  file sample
  

  //the --id is = to --sample. this way the dir that is being created is the same as the sample name for each sample.
  //the localmem and cores specifies how much memory and cpu the task can use
  """
  cellranger count \
  --id=${sample} \
  --transcriptome=${tx_reference} \
  --fastqs=${fastq_dir} \
  --sample=${sample} \
  --localcores ${task.cpus} \
  --localmem ${task.memory.toGiga()}

  """
}