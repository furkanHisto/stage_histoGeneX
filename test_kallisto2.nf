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
        .groupTuple(by:[0])                                                             // you then sort by tuple. you can specify by which tuple value you want to sort. because we specified in the earlier step that sample\d is the first one, it will sort by sample by default.
        .set { read_files_kallisto}
        

process prep {
        input:
        set val(name), file(reads) from read_files_kallisto
        output:
        set val(name), file(reads)) into preprocess

        
        """
        cat  ${reads}
        echo ${reads}
        """
}

preprocess.view()

// process kallisto_bus {
//     tag "${name}"

//     input: 
// //     set val(name), file(reads) from read_files_kallisto
//     tuple val(name), file (reads) from preprocess
//     file index from kallisto_index.collect()


//     """
//     kallisto bus  -i $index \
//                     --output=/home/histogenex/test_kallisto/test1/${name} \
//                     -x '10xv3' -t ${task.cpus} \
//                     ${reads} | tee ${name}_kallisto.log
//     """
// }
