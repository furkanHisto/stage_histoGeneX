#!/usr/bin/env nextflow

transcriptome_reference = Channel
                                .fromPath(params.transcriptome)
                                .ifEmpty { exit 1, "Reference not found: ${params.transcriptome}" }



fastqdir_channel = Channel
                          .fromPath(params.fastq_dir)
                          .ifEmpty { exit 1, "Samples not found: ${params.fasq_dir}" }


sampleChannel = Channel.fromPath("${params.fastq_dir}/*.gz")
  .map { Path path -> 
    path.toFile()
        .getSimpleName()                                       // this returns the name without the extension
        .replaceAll(/_S[1-4]_L00[1-2]_[IR][1-2]_001/, '')      // this replaces the first string with something else. we need to change a part of the file name so the cellranger count can recognize the file name
  }
  .unique()  //this makes it so that there are no duplicate files



process CellRangerCount {

  publishDir '${params.outdir}/cellrangercount',  mode: 'copy'       // this is added to specify the directory of the output of the run.

  //label "mid_memory"     // this looks in the hgx.config file for this label. this label sets the configuration of all the processes to set values.
  tag "$sample"         // this is added to check wich sample is currently running when executing nextflow.

  input:
  val sample from sampleChannel      // we take the values generated from the samplechannel. these are the file names 'GEM_sample*' 
  file fastq_dir from fastqdir_channel.first()              // we add .first because that channel only outputs 1 file while the samplechannel outputs 4. when there is no .first() then the run will stop at 1 file because the fastq_dir channel only outputs 1 file. 
  file tx_reference from transcriptome_reference.first()    // .first takes the first output file and reuses it. so the only output of these channels is reuses 4 times to match the amount of samplechannel output.
  
  

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