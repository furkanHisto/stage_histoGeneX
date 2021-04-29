#!/usr/bin/env nextflow

kallisto_index = Channel.fromPath(params.kallisto_index)
                        .ifEmpty { exit 1, "Index not found: ${params.kallisto_index}" }

fastqdir_channel = Channel
                                .fromPath(params.fastq_dir)
                                .ifEmpty { exit 1, "Samples not found: ${params.fasq_dir}" }
                            
kallisto_sampleChannel = Channel.fromPath("${params.fastq_dir}/*.gz")
    .reduce([]){                    //it accumulates the files and makes 1 list. first you make an empty array to put the values in [].
        acc,item ->                 // define acc and item
            acc << item             // all the values of the channel are accumulated. then you return the value into the empty array.
            return acc
            
    }
    .map{
        list -> list.join(" ")      // you join the whole array/list in 1 string that is space seperated. this way you can put all the files into the kallisto command.
    }

output_kallisto = Channel.fromPath("${params.output_kallisto}")

process kallisto_quant {
    tag "${sample}"

    input: 
    val sample from kallisto_sampleChannel

    """
    kallisto quant  --index=${params.kallisto_index} \
                    --output=${params.output_kallisto} \
                    ${sample} 
    """
}