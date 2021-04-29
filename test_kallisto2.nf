#!/usr/bin/env nextflow

//creating channel to get the kallisto index file
kallisto_index = Channel.fromPath(params.kallisto_index)
                        .ifEmpty { exit 1, "Index not found: ${params.kallisto_index}" }

// this method is not good. that's why it's in comment                            
// kallisto_sampleChannel = Channel.fromPath("${params.fastq_dir}/*.gz")
//                                .map{                                
//                                     Path path ->
//                                     path.toFile()
//                                     .getAbsolutePath()              // wih this the absolute path of the file is given. because kallisto quant needs the paths of the file, we add this
//                                } 
//                                .collate(2)                          // with collate the emitted values are grouped in tuples. the 2 indicates how many values are combined into a tuple. 


//channel to read the sample files
Channel.fromFilePairs(params.input).view()
        .ifEmpty {exit 1, "params.input_paths was empty - no input files supplied" }
        .map {tag, pair -> 
                subtags = (tag =~ /sample\d{1}/) [0];
                tuple ( subtags, pair )
        }
        
        .groupTuple(by:0).view()

        // .map {tag, pair -> 
        //         subtags = (tag =~ /L00\d{1}/) [0];
        //         tuple [subtags[1], pair]
        // }
        // .groupTuple(by:0).view()

        // .set { read_files_kallisto}


             

/*
process kallisto_bus {
    tag "${name}"

    input: 
    set val(name), file(reads) from read_files_kallisto
    file index from kallisto_index.collect()

    """
    kallisto bus  -i $index \
                    --output=/home/histogenex/test_kallisto/test1/${name} \
                    -x '10xv3' -t ${task.cpus} \
                    $reads | tee ${name}_kallisto.log
    """
}
*/