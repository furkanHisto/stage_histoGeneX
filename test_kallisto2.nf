#!/usr/bin/env nextflow


//creating channel to get the kallisto index file
kallisto_index = Channel.fromPath(params.kallisto_index)
                        .ifEmpty { exit 1, "Index not found: ${params.kallisto_index}" }



//channel to read the sample files
Channel.fromFilePairs(params.input)                                                     // the files emited are collected in pairs into a tuple. 
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

        
        



process kallisto{
        tag "${name}"
        input:
        tuple val (name), file (read1), file (read2) from read_files_kallisto
        file index from kallisto_index.collect()
        output:


         """
         kallisto bus  -i $index \
                       --output=/home/histogenex/test_kallisto/test2/${name} \
                       -x '10xv2' -t ${task.cpus} \
                       ${read1} ${read2} | tee ${name}_kallisto.log
        """
}



