#!/usr/bin/env nextflow


sampleChannel = Channel.fromPath('/mnt/storage/Data/10xGenomics_TrainingsRun/*.gz')
  .map { Path path -> 
    path.toFile()
        .replaceAll(/_S[1-4]_L00[1-2]_[IR][1-2]_001.fastq.gz/, '')      // this replaces the first string with something else. we need to change a part of the file name so the cellranger count can recognize the file name
  }
  .unique()  //this makes it so that every duplicate file stays as 1 unique file

process CellRangeCount {
  input:
  val sample from sampleChannel

  """
  cellranger count --id=cellrangerCount --transcriptome=/mnt/storage/Reference/refdata-gex-GRCh38-2020-A --fastqs=/mnt/storage/Data/10xGenomics_TrainingsRun --sample=${sample}
  """
}