/*
 * COMPLETE PIPELINE: 
 * Created by Furkan
 * Changes by HGX, PJ
 */


/*
 * USAGE
 */


// The typical command for running the pipeline is as follows:

// nextflow run hgx_rnaseq.nf -profile hgx,docker

// Mandatory arguments:
// -profile   Configuration profile to use. Can use multiple (comma seperated)


/*
 * Set up Configuration Variables
 */


// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the Genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}


// Reference index path configuration
// Define these here - after the profiles are loaded with the Genomes paths (from genomes.config file)
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.transcript_fasta = params.genome ? params.genomes[ params.genome ].transcript_fasta ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.gtf_ref_flat = params.genome ? params.genomes[ params.genome ].gtf_ref_flat ?: false : false
params.tx2gene = params.genome ? params.genomes[ params.genome ].tx2gene ?: false : false
params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
params.kallisto_index = params.genome ? params.genomes[ params.genome ].kallisto_index ?: false : false
params.star_index = params.genome ? params.genomes[ params.genome ].star_index ?: false : false
params.rRNA_database_manifest = params.genome ? params.genomes[ params.genome ].rRNA_database_manifest ?: false : false
params.blacklist = params.genome ? params.genomes[ params.genome ].blacklist ?: false : false
params.cytoband = params.genome ? params.genomes[ params.genome ].cytoband ?: false : false
params.protDom = params.genome ? params.genomes[ params.genome ].protDom ?: false : false
params.multiqc_config = params.genome ? params.genomes[ params.genome ].multiqc_config ?: false : false


ch_mdsplot_header = Channel.fromPath("$baseDir/assets/mdsplot_header.txt", checkIfExists: true)
ch_heatmap_header = Channel.fromPath("$baseDir/assets/heatmap_header.txt", checkIfExists: true)
ch_biotypes_header = Channel.fromPath("$baseDir/assets/biotypes_header.txt", checkIfExists: true)
Channel.fromPath("$baseDir/assets/where_are_my_files.txt", checkIfExists: true)
       .into{ch_where_trim_galore; ch_where_star}


// Define regular variables so that they can be overwritten
clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_clip_r2 = params.three_prime_clip_r2


// Get rRNA databases
rRNA_database = file(params.rRNA_database_manifest)
if (rRNA_database.isEmpty()) {exit 1, "File ${rRNA_database.getName()} is empty!"}
Channel
    .from( rRNA_database.readLines() )
    .map { row -> file(row) }
    .set { sortmerna_fasta }
  

// Validate inputs
if (params.star_index && params.aligner == 'star' && !params.skipAlignment) {
  if (hasExtension(params.star_index, 'gz')) {
    star_index_gz = Channel
        .fromPath(params.star_index, checkIfExists: true)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
  } else{
    star_index = Channel
        .fromPath(params.star_index, checkIfExists: true)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
  }
}


if (params.fasta) {
  fasta = Channel
      .fromPath(params.fasta, checkIfExists: true)
      .ifEmpty { exit 1, "FASTA file not found: ${params.fasta}" }
}


if (params.blacklist) {
  blacklist = Channel
      .fromPath(params.blacklist, checkIfExists: true)
      .ifEmpty { exit 1, "Blacklist not found: ${params.blacklist}" }
}


if (params.cytoband) {
  cytoband = Channel
      .fromPath(params.cytoband, checkIfExists: true)
      .ifEmpty { exit 1, "Cytoband not found: ${params.cytoband}" }
}


if (params.protDom) {
  protDom = Channel
      .fromPath(params.protDom, checkIfExists: true)
      .ifEmpty { exit 1, "ProtDom not found: ${params.protDom}" }
}


if (params.gtf_ref_flat) {
  gtf_ref_flat = Channel
      .fromPath(params.gtf_ref_flat, checkIfExists: true)
      .ifEmpty { exit 1, "GTF REF_FLAT file not found: ${params.gtf_ref_flat}" }
}


if (params.tx2gene) {
  tx2gene_tximport = Channel
      .fromPath(params.tx2gene, checkIfExists: true)
      .ifEmpty { exit 1, "GENCODE CSV file not found: ${params.tx2gene}" }
}


if (params.tx2gene) {
  tx2gene_merge = Channel
      .fromPath(params.tx2gene, checkIfExists: true)
      .ifEmpty { exit 1, "GENCODE CSV file not found: ${params.tx2gene}" }
}


Channel
  .fromPath(params.fasta, checkIfExists: true)
    .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
    .set { genome_fasta_gz }


// Check whether kallisto needs a transcriptome fasta to exctract transcripts from.
if (params.pseudo_aligner == 'kallisto') {
  if (params.kallisto_index) {
    if (hasExtension(params.kallisto_index, 'gz')) {
      kallisto_index_gz = Channel
          .fromPath(params.kallisto_index, checkIfExists: true)
          .ifEmpty { exit 1, "kallisto index not found: ${params.kallisto_index}" }
    } else {
      kallisto_index = Channel
          .fromPath(params.kallisto_index, checkIfExists: true)
          .ifEmpty { exit 1, "kallisto index not found: ${params.kallisto_index}" }
      }
  }
}

// Set transcriptome fasta
if (hasExtension(params.transcript_fasta, 'gz')) {
  transcript_fasta_gz = Channel
      .fromPath(params.transcript_fasta, checkIfExists: true)
      .ifEmpty { exit 1, "Transcript fasta file not found: ${params.transcript_fasta}" }
} else {
  ch_fasta_for_kallisto_index = Channel
      .fromPath(params.transcript_fasta, checkIfExists: true)
      .ifEmpty { exit 1, "Transcript fasta file not found: ${params.transcript_fasta}" }
  }


if (params.gtf) {
  if (params.gff) {
      // Prefer gtf over gff
      log.info "Both GTF and GFF have been provided: Using GTF as priority."
  }
  if (hasExtension(params.gtf, 'gz')) {
  gtf_gz = Channel
        .fromPath(params.gtf, checkIfExists: true)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
  } else {
    Channel
        .fromPath(params.gtf, checkIfExists: true)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtf_makeSTARindex; gtf_makekallistoIndex; gtf_makeBED12;
                gtf_star; gtf_dupradar; gtf_featureCounts; gtf_kallisto; gtf_kallisto_merge; gtf_arriba; gtf_fusions }

    }
} else if (params.gff) {
  if (hasExtension(params.gff, 'gz')) {
    gff_gz = Channel.fromPath(params.gff, checkIfExists: true)
                  .ifEmpty { exit 1, "GFF annotation file not found: ${params.gff}" }
  } else {
    gffFile = Channel.fromPath(params.gff, checkIfExists: true)
                  .ifEmpty { exit 1, "GFF annotation file not found: ${params.gff}" }
    }
  } else {
    exit 1, "No GTF or GFF3 annotation specified!"
    }


if (params.bed12) {
    bed12 = Channel
        .fromPath(params.bed12, checkIfExists: true)
        .ifEmpty { exit 1, "BED12 annotation file not found: ${params.bed12}" }
        .set { bed_rseqc }
}

// Gencode parameter for feature counts biotypes
if (params.gencode) {
  biotype = "gene_type"
} else {
  biotype = params.fc_group_features_type
  }


// Stage config files
ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)


// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}


/*
 * Create a channel for input read files
 */


if (params.readPaths) {
  if (params.singleEnd) {
      Channel
          .from(params.readPaths)
          .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
          .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
          .into { raw_reads_fastqc; raw_reads_trimgalore }
  } else {
      Channel
          .from(params.readPaths)
          .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ] }
          .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
          .into { raw_reads_fastqc; raw_reads_trimgalore }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { raw_reads_fastqc; raw_reads_trimgalore }
  }


// Unzip fasta and gtf/gff files
process gunzip_genome_fasta {
        tag "$gz"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gz from genome_fasta_gz

        output:
        file "${gz.baseName}" into ch_fasta_for_star_index, ch_fasta_for_kallisto_transcripts

        script:
        """
        gunzip --verbose --stdout --force ${gz} > ${gz.baseName}
        """
}


process gunzip_gtf {
        tag "$gz"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gz from gtf_gz

        output:
        file "${gz.baseName}" into gtf_makeSTARindex, gtf_makekallistoIndex, gtf_makeBED12,
                                        gtf_star, gtf_dupradar, gtf_featureCounts, gtf_kallisto, gtf_kallisto_merge, gtf_arriba, gtf_fusions

        script:
        """
        gunzip --verbose --stdout --force ${gz} > ${gz.baseName}
        """
}


if (params.gff && !params.gtf) {
    process gunzip_gff {
        tag "$gz"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                    saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gz from gff_gz

        output:
        file "${gz.baseName}" into gffFile

        script:
        """
        gunzip --verbose --stdout --force ${gz} > ${gz.baseName}
        """
    }
}


if (params.transcript_fasta && params.pseudo_aligner == 'kallisto' && !params.kallisto_index) {
    process gunzip_transcript_fasta {
        tag "$gz"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_transcriptome" : params.outdir },
                    saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gz from transcript_fasta_gz

        output:
        file "${gz.baseName}" into ch_fasta_for_kallisto_index

        script:
        """
        gunzip --verbose --stdout --force ${gz} > ${gz.baseName}
        """
    }
}


/*
 * PREPROCESSING - Build BED12 file
 */

if (!params.bed12) {
    process makeBED12 {
        tag "$gtf"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                    saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gtf from gtf_makeBED12

        output:
        file "${gtf.baseName}.bed" into bed_rseqc

        script: // This script is bundled with the pipeline, inserted in docker.
        """
        gtf2bed $gtf > ${gtf.baseName}.bed
        """
    }
}


/*
 * PREPROCESSING - Build STAR index
 */

if (!params.skipAlignment) {
  if (params.aligner == 'star' && !params.star_index && params.fasta) {
      process makeSTARindex {
          label 'high_memory'
          tag "$fasta"
          publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                      saveAs: { params.saveReference ? it : null }, mode: 'copy'

          input:
          file fasta from ch_fasta_for_star_index
          file gtf from gtf_makeSTARindex

          output:
          file 'star' into star_index

          script:
          def avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
          """
          mkdir star
          STAR \\
              --runMode genomeGenerate \\
              --runThreadN 10 \\
              --sjdbGTFfile $gtf \\
              --genomeDir star/ \\
              --genomeFastaFiles $fasta \\
              $avail_mem
          """
      }
  } 
}


/*
 * PREPROCESSING - Create Kallisto transcriptome index
 */

if (params.pseudo_aligner == 'kallisto' && !params.kallisto_index) {
    process makekallistoIndex {
        label "kallisto"
        tag "$fasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                            saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from ch_fasta_for_kallisto_index

        output:
        file 'kallisto_index' into kallisto_index

        script:
        //def gencode = params.gencode  ? '--gencode' : ''
        """
        kallisto index $fasta -i kallisto_index
        """
    }
}


/*
 * STEP 1 - FastQC
 */
 
 // Quality Control - Checks on raw sequence data.


 process fastqc {
    tag "$name"
    label 'low_memory'
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename" }

    when:
    !params.skipQC && !params.skipFastQC

    input:
    set val(name), file(reads) from raw_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc --quiet --threads 2 $reads
    """
}


/*
 * STEP 2 - Trim Galore!
 */

// Quality and Adapter trimming with FASTQC on trimmed FASTQ.


if (!params.skipTrimming) {
    process trim_galore {
        label 'low_memory'
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else if (!params.saveTrimmed && filename == "where_are_my_files.txt") filename
                else if (params.saveTrimmed && filename != "where_are_my_files.txt") filename
                else null
            }

        input:
        set val(name), file(reads) from raw_reads_trimgalore
        file wherearemyfiles from ch_where_trim_galore.collect()

        output:
        set val(name), file("*fq.gz") into trimgalore_reads
        file "*trimming_report.txt" into trimgalore_results
        file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
        file "where_are_my_files.txt"

        script:
        c_r1 = clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
        c_r2 = clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
        tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
        tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''
        nextseq = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ''
        if (params.singleEnd) {
            """
            trim_galore --fastqc --gzip $c_r1 $tpc_r1 $nextseq $reads
            """
        } else {
            """
            trim_galore --paired --fastqc --gzip -j 2 $c_r1 $c_r2 $tpc_r1 $tpc_r2 $nextseq $reads
            """
          }
    }
} else {
   raw_reads_trimgalore
       .set {trimgalore_reads}
   trimgalore_results = Channel.empty()
  }


/*
 * STEP 2+ - SortMeRNA 
 */
 
// Remove rRNA sequences on request.


if (!params.removeRiboRNA) {
    trimgalore_reads
        .into { trimmed_reads_alignment; trimmed_reads_kallisto }
    sortmerna_logs = Channel.empty()
} else {
    process sortmerna_index {
        label 'low_memory'
        tag "${fasta.baseName}"

        input:
        file(fasta) from sortmerna_fasta

        output:
        val("${fasta.baseName}") into sortmerna_db_name
        file("$fasta") into sortmerna_db_fasta
        file("${fasta.baseName}*") into sortmerna_db

        script:
        """
        indexdb_rna --ref $fasta,${fasta.baseName} -m 3072 -v
        """
    }

    process sortmerna {
        label 'low_memory'
        tag "$name"
        publishDir "${params.outdir}/SortMeRNA", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_rRNA_report.txt") > 0) "logs/$filename"
                else if (params.saveNonRiboRNAReads) "reads/$filename"
                else null
            }

        input:
        set val(name), file(reads) from trimgalore_reads
        val(db_name) from sortmerna_db_name.collect()
        file(db_fasta) from sortmerna_db_fasta.collect()
        file(db) from sortmerna_db.collect()

        output:
        set val(name), file("*.fq.gz") into trimmed_reads_alignment, trimmed_reads_kallisto
        file "*_rRNA_report.txt" into sortmerna_logs


        script:
        //concatenate reference files: ${db_fasta},${db_name}:${db_fasta},${db_name}:...
        def Refs = ''
        for (i=0; i<db_fasta.size(); i++) { Refs+= ":${db_fasta[i]},${db_name[i]}" }
        Refs = Refs.substring(1)

        if (params.singleEnd) {
            """
            gzip -d --force < ${reads} > all-reads.fastq

            sortmerna --ref ${Refs} \
                --reads all-reads.fastq \
                --num_alignments 1 \
                -a ${task.cpus} \
                --fastx \
                --aligned rRNA-reads \
                --other non-rRNA-reads \
                --log -v

            gzip --force < non-rRNA-reads.fastq > ${name}.fq.gz

            mv rRNA-reads.log ${name}_rRNA_report.txt
            """
        } else {
            """
            gzip -d --force < ${reads[0]} > reads-fw.fq
            gzip -d --force < ${reads[1]} > reads-rv.fq
            merge-paired-reads.sh reads-fw.fq reads-rv.fq all-reads.fastq

            sortmerna --ref ${Refs} \
                --reads all-reads.fastq \
                --num_alignments 1 \
                -a ${task.cpus} \
                --fastx --paired_in \
                --aligned rRNA-reads \
                --other non-rRNA-reads \
                --log -v

            unmerge-paired-reads.sh non-rRNA-reads.fastq non-rRNA-reads-fw.fq non-rRNA-reads-rv.fq
            gzip < non-rRNA-reads-fw.fq > ${name}-fw.fq.gz
            gzip < non-rRNA-reads-rv.fq > ${name}-rv.fq.gz

            mv rRNA-reads.log ${name}_rRNA_report.txt
            """
          }
    }
}  


/*
 * STEP 3 - align with STAR
 */
 
// STAR is an ultrafast universal RNA-seq aligner.



// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false.

skipped_poor_alignment = []
def check_log(logs) {
    def percent_aligned = 0;
    logs.eachLine { line ->
        if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
            percent_aligned = matcher[0][1]
        }
    }
    logname = logs.getBaseName() - 'Log.final'
    if (percent_aligned.toFloat() <= '5'.toFloat()) {
        log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percent_aligned}% <<"
        skipped_poor_alignment << logname
        return false
    } else {
        log.info "          Passed alignment > star ($logname)   >> ${percent_aligned}% <<"
        return true
      }
}



// STAR function
// Small adaptions for arriba


if (!params.skipAlignment) {
  if (params.aligner == 'star') {
      hisat_stdout = Channel.from(false)
      process star {
          label 'high_memory'
          tag "$name"
          publishDir "${params.outdir}/STAR", mode: 'copy',
              saveAs: {filename ->
                  if (filename.indexOf(".bam") == -1) "logs/$filename"
                  else if (params.saveUnaligned && filename != "where_are_my_files.txt" && 'Unmapped' in filename) unmapped/filename
                  else if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
                  else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") filename
                  else null
              }

          input:
          set val(name), file(reads) from trimmed_reads_alignment
          file index from star_index.collect()
          file gtf from gtf_star.collect()
          file wherearemyfiles from ch_where_star.collect()

          output:
          set file("*Log.final.out"), file ('*.bam') into star_aligned
          file "*.out" into alignment_logs
          file "*SJ.out.tab"
          file "*Log.out" into star_log
          file "where_are_my_files.txt"
          file "*Unmapped*" optional true
          file "${prefix}Aligned.sortedByCoord.out.bam.bai" into bam_index_rseqc, bam_index_genebody, bai_fusions

          script:
          prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
          def star_mem = task.memory ?: params.star_memory ?: false
          def avail_mem = star_mem ? "--limitBAMsortRAM ${star_mem.toBytes() - 100000000}" : ''
          seq_center = params.seq_center ? "--outSAMattrRGline ID:$prefix 'CN:$params.seq_center' 'SM:$prefix'" : "--outSAMattrRGline ID:$prefix 'SM:$prefix'"
          unaligned = params.saveUnaligned ? "--outReadsUnmapped Fastx" : ''
          """
          STAR --genomeDir $index \\
              --sjdbGTFfile $gtf \\
              --readFilesIn $reads  \\
              --runThreadN ${task.cpus} \\
              --twopassMode Basic \\
              --outWigType bedGraph \\
              --outSAMtype BAM SortedByCoordinate $avail_mem \\
              --readFilesCommand zcat \\
              --runDirPerm All_RWX $unaligned \\
              --outFileNamePrefix $prefix $seq_center \\
			  --chimSegmentMin $params.ChimSegMin \\
			  --chimOutType WithinBAM SoftClip \\
			  --chimJunctionOverhangMin $params.OverhangMin \\
			  --chimScoreMin $params.ScoreMin \\
			  --chimScoreDropMax $params.ScoreDropMax \\
			  --chimScoreJunctionNonGTAG $params.JunctionNonGTAG \\
			  --chimScoreSeparation $params.ScoreSeparation \\
			  --alignSJstitchMismatchNmax $params.MismatchNmax \\
			  --chimSegmentReadGapMax $params.ReadGapMax

          samtools index ${prefix}Aligned.sortedByCoord.out.bam
          """
      }

      // Filter removes all 'aligned' channels that fail the check.

      star_aligned
          .filter { logs, bams -> check_log(logs) }
          .flatMap {  logs, bams -> bams }
      .into { bam_count; bam_rseqc; bam_preseq; bam_markduplicates; bam_featurecounts; bam_arriba; bam_rnaseqmetrics; bam_fusions }
  }
}


/*
 * STEP 4 - Arriba
 */

// Fast and accurate gene fusion detection from RNA-seq data.


if (!params.SkipArriba) {
    process arriba {
      label 'low_memory'
      tag "${bam.simpleName}"
      publishDir "${params.outdir}/Arriba", mode: 'copy'

      input:
      file bam from bam_arriba
      file gtf from gtf_arriba.collect()
      file fasta from fasta.collect()
      file blacklist from blacklist.collect()
	  
      output:
      file "${bam.baseName}.fusions.tsv" into tsvfusions
      file "${bam.baseName}.fusions.discarded.tsv" into tsvfiles
	  
      script:
      """
      arriba -x $bam \\
        -o ${bam.baseName}.fusions.tsv \\
		    -O ${bam.baseName}.fusions.discarded.tsv \\
		    -a $fasta \\
		    -g $gtf \\
		    -b $blacklist \\
		    -T -P
      """
    }
}

// Fusionvisualizer to create pdfs of fusions
// Incorporated in arriba

if (!params.SkipArriba) {
    process FusionVisualizer {
      label 'low_memory'
      tag "${fusionbam.simpleName}"
      publishDir "${params.outdir}/Arriba", mode: 'copy'

      input:
      file fusionbam from bam_fusions
      file gtf from gtf_fusions.collect()
      file tsv from tsvfusions
      file bai from bai_fusions
      file cytoband from cytoband.collect()
      file protDom from protDom.collect()

      output:
      file "${fusionbam.baseName}.fusions.pdf" into fusionspdf
	  
      script:
      """
      draw_fusions.R --annotation=$gtf \\
		    --fusions=$tsv \\
		    --output=${fusionbam.baseName}.fusions.pdf \\
		    --alignments=$fusionbam \\
		    --cytobands=$cytoband \\
        --proteinDomains=$protDom
      """
    }
}


/*
 * STEP 5 - RSeQC analysis
 */

// RNA-seq quality control, provides a number of useful modules that can comprehensively evaluate high throughput RNA-seq data.


  process rseqc {
      label 'mid_memory'
      tag "${bam_rseqc.baseName - '.sorted'}"
      publishDir "${params.outdir}/rseqc" , mode: 'copy',
          saveAs: {filename ->
                   if (filename.indexOf("bam_stat.txt") > 0)                      "bam_stat/$filename"
              else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
              else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
              else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
              else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
              else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
              else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
              // else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
              // else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
              // else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
              // else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
              else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
              else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
              else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
              else if (filename.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$filename"
              else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
              else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
              else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
              else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
              else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
              else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
              else filename
          }

      when:
      !params.skipQC && !params.skipRseQC

      input:
      file bam_rseqc
      file index from bam_index_rseqc
      file bed12 from bed_rseqc.collect()

      output:
      file "*.{txt,pdf,r,xls}" into rseqc_results

      script:
      """
      infer_experiment.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.infer_experiment.txt
      junction_annotation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12 2> ${bam_rseqc.baseName}.junction_annotation_log.txt
      bam_stat.py -i $bam_rseqc > ${bam_rseqc.baseName}.bam_stat.txt
      junction_saturation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12 
      inner_distance.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
      read_distribution.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.read_distribution.txt
      read_duplication.py -i $bam_rseqc -o ${bam_rseqc.baseName}.read_duplication
      """
      // RPKM_saturation.py -r $bed12 -d "1+-,1-+,2++,2--" -i $bam_rseqc -o ${bam_rseqc.baseName}.RPKM_saturation

  }
  
    
/*
 * STEP 6 - preseq analysis
 */
  
// Preseq estimates the complexity of a library.


  process preseq {
      tag "${bam_preseq.baseName - '.sorted'}"
      publishDir "${params.outdir}/preseq", mode: 'copy'

      when:
      !params.skipQC && !params.skipPreseq

      input:
      file bam_preseq

      output:
      file "${bam_preseq.baseName}.ccurve.txt" into preseq_results

      script:
      """
      preseq lc_extrap -v -B $bam_preseq -o ${bam_preseq.baseName}.ccurve.txt
      """
  }


/*
 * STEP 7 - Mark duplicates
 */

// Identifies duplicate reads.


  process markDuplicates {
      tag "${bam.baseName - '.sorted'}"
      publishDir "${params.outdir}/markDuplicates", mode: 'copy',
          saveAs: {filename -> filename.indexOf("_metrics.txt") > 0 ? "metrics/$filename" : "$filename"}

      when:
      !params.skipQC && !params.skipDupRadar

      input:
      file bam from bam_markduplicates

      output:
      file "${bam.baseName}.markDups.bam" into bam_md
      file "${bam.baseName}.markDups_metrics.txt" into picard_results
      file "${bam.baseName}.markDups.bam.bai"

      script:
      markdup_java_options = (task.memory.toGiga() > 8) ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2 )+"g "+ "-Xmx" + (task.memory.toGiga() - 1)+ "g\""
      """
      picard ${markdup_java_options} MarkDuplicates \\
          INPUT=$bam \\
          OUTPUT=${bam.baseName}.markDups.bam \\
          METRICS_FILE=${bam.baseName}.markDups_metrics.txt \\
          REMOVE_DUPLICATES=false \\
          ASSUME_SORTED=true \\
          PROGRAM_RECORD_ID='null' \\
          VALIDATION_STRINGENCY=LENIENT
      samtools index ${bam.baseName}.markDups.bam
      """
  }


/*
 * STEP 8 - RnaSeqMetrics
 */

// RnaSeqMetrics calculates standard RNA-seq related metrics.


   process RnaSeqMetrics {
		tag "${bam.baseName - '.sorted'}"
		publishDir "${params.outdir}/RnaSeqMetrics", mode: 'copy'
	
		when:
		!params.skipQC && !params.skipRnaSeqMetrics
	
		input:
		file bam from bam_rnaseqmetrics
		file ref_flat from gtf_ref_flat.collect()
	
		output:
		file "${bam.baseName}.rna_metrics" into rnaseqmetrics_results
	
		script:
		rnaseqmetrics_java_options = (task.memory.toGiga() > 8) ? params.rnaseqmetrics_java_options : "\"-Xms" +  (task.memory.toGiga() / 2 )+"g "+ "-Xmx" + (task.memory.toGiga() - 1)+ "g\""
		"""
		picard ${rnaseqmetrics_java_options} CollectRnaSeqMetrics \\
			INPUT=$bam \\
			OUTPUT=${bam.baseName}.rna_metrics \\
			REF_FLAT=$ref_flat \\
			STRAND=SECOND_READ_TRANSCRIPTION_STRAND \\
			RIBOSOMAL_INTERVALS='null'
		"""
   }


/*
 * STEP 9 - dupRadar
 */

// DupRadar provides duplication rate quality control for RNA-seq datasets.


  process dupradar {
      label 'low_memory'
      tag "${bam_md.baseName - '.sorted.markDups'}"
      publishDir "${params.outdir}/dupradar", mode: 'copy',
          saveAs: {filename ->
              if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
              else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
              else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
              else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
              else if (filename.indexOf("_duprateExpDensCurve.txt") > 0) "scatter_curve_data/$filename"
              else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
              else "$filename"
          }

      when:
      !params.skipQC && !params.skipDupRadar

      input:
      file bam_md
      file gtf from gtf_dupradar.collect()

      output:
      file "*.{pdf,txt}" into dupradar_results

      script: // This script is bundled with the pipeline, inserted in docker.
      def dupradar_direction = 0
      if (params.forwardStranded && !params.unStranded) {
          dupradar_direction = 1
      } else if (params.reverseStranded && !params.unStranded) {
          dupradar_direction = 2
      }
      def paired = params.singleEnd ? 'single' :  'paired'
      """
      dupRadar.r $bam_md $gtf $dupradar_direction $paired ${task.cpus}
      """
  }


/*
 * STEP 10a - Feature counts
 */

 // Counts mapped reads for genomic features such as genes, exons, promotor, gene bodies, genomic bins and chromosomal locations.
  

  process featureCounts {
      label 'low_memory'
      tag "${bam_featurecounts.baseName - '.sorted'}"
      publishDir "${params.outdir}/featureCounts", mode: 'copy',
          saveAs: {filename ->
              if (filename.indexOf("biotype_counts") > 0) "biotype_counts/$filename"
              else if (filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
              else if (filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
              else "$filename"
          }

      input:
      file bam_featurecounts
      file gtf from gtf_featureCounts.collect()
      file biotypes_header from ch_biotypes_header.collect()

      output:
      file "${bam_featurecounts.baseName}_gene.featureCounts.txt" into geneCounts, featureCounts_to_merge
      file "${bam_featurecounts.baseName}_gene.featureCounts.txt.summary" into featureCounts_logs
      file "${bam_featurecounts.baseName}_biotype_counts*mqc.{txt,tsv}" optional true into featureCounts_biotype

      script:
      def featureCounts_direction = 0
      def extraAttributes = params.fc_extra_attributes ? "--extraAttributes ${params.fc_extra_attributes}" : ''
      if (params.forwardStranded && !params.unStranded) {
          featureCounts_direction = 1
      } else if (params.reverseStranded && !params.unStranded) {
          featureCounts_direction = 2
      }
      // Try to get real sample name
      sample_name = bam_featurecounts.baseName - '_R1_001' - 'Aligned.sortedByCoord.out' - '_subsamp.sorted' 
      biotype_qc = params.skipBiotypeQC ? '' : "featureCounts -a $gtf -g $biotype -o ${bam_featurecounts.baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts"
      mod_biotype = params.skipBiotypeQC ? '' : "cut -f 1,7 ${bam_featurecounts.baseName}_biotype.featureCounts.txt | tail -n +3 | cat $biotypes_header - >> ${bam_featurecounts.baseName}_biotype_counts_mqc.txt && mqc_features_stat.py ${bam_featurecounts.baseName}_biotype_counts_mqc.txt -s $sample_name -f rRNA -o ${bam_featurecounts.baseName}_biotype_counts_gs_mqc.tsv"
      """
      featureCounts -T ${task.cpus} -a $gtf -g ${params.fc_group_features} -t ${params.fc_count_type} -o ${bam_featurecounts.baseName}_gene.featureCounts.txt $extraAttributes -p -s $featureCounts_direction $bam_featurecounts
      $biotype_qc
      $mod_biotype
      """
  }


/*
 * STEP 10b - Merge featurecounts for all samples
 */


  process merge_featureCounts {
      label 'high_memory'
      tag "${input_files[0].baseName - '.sorted'}"
      publishDir "${params.outdir}/featureCounts", mode: 'copy'

      input:
      file input_files from featureCounts_to_merge.collect()

      output:
      file 'merged_gene_counts.txt' into featurecounts_merged

      script:
      // Redirection (the `<()`) for the win!
      // Geneid in 1st column and gene_name in 7th
      gene_ids = "<(tail -n +2 ${input_files[0]} | cut -f1,7 )"
      counts = input_files.collect{filename ->
        // Remove first line and take third column
        "<(tail -n +2 ${filename} | sed 's:.sorted.bam::' | cut -f8)"}.join(" ")
      """
      paste $gene_ids $counts > merged_gene_counts.txt
      """
  }


/*
 * STEP 11 - edgeR MDS and heatmap
 */

// Sample similarity is generated from mormalised gene counts through edgeR.
  

  process sample_correlation {
      label 'low_memory'
      tag "${input_files[0].toString() - '.sorted_gene.featureCounts.txt' - 'Aligned'}"
      publishDir "${params.outdir}/sample_correlation", mode: 'copy'

      when:
      !params.skipQC && !params.skipEdgeR

      input:
      file input_files from geneCounts.collect()
      val num_bams from bam_count.count()
      file mdsplot_header from ch_mdsplot_header.collect()
      file heatmap_header from ch_heatmap_header.collect()

      output:
      file "*.{txt,pdf,csv}" into sample_correlation_results

      when:
      num_bams > 2 && (!params.sampleLevel)

      script: // This script is bundled with the pipeline, inserted in docker.
      """
      edgeR_heatmap_MDS.r $input_files
      cat $mdsplot_header edgeR_MDS_Aplot_coordinates_mqc.csv >> tmp_file
      mv tmp_file edgeR_MDS_Aplot_coordinates_mqc.csv
      cat $heatmap_header log2CPM_sample_correlation_mqc.csv >> tmp_file
      mv tmp_file log2CPM_sample_correlation_mqc.csv
      """
  }


/*
 * STEP 12a - Transcripts quantification with kallisto
 */

// Kallisto is a program for quantifying abundances of target sequences using high-throughput sequencing reads.


 if (params.pseudo_aligner == 'kallisto') {
    process kallisto {
        label 'kallisto'
        tag "$sample"
        publishDir "${params.outdir}/kallisto", mode: 'copy'

        input:
        set sample, file(reads) from trimmed_reads_kallisto
        file index from kallisto_index.collect()
        file gtf from gtf_kallisto.collect()

        output:
        file "${sample}/" into kallisto_logs
        set val(sample), file("${sample}/") into kallisto_tximport, kallisto_parsegtf

        script:
        def rnastrandness = params.singleEnd ? 'U' : 'IU'
        if (params.forwardStranded && !params.unStranded) {
            rnastrandness = params.singleEnd ? 'SF' : 'ISF'
        } else if (params.reverseStranded && !params.unStranded) {
            rnastrandness = params.singleEnd ? 'SR' : 'ISR'
        }
        def endedness = params.singleEnd ? "-r ${reads[0]}" : "${reads[0]} ${reads[1]}"
        unmapped = params.saveUnaligned ? "--writeUnmappedNames" : ''
        """
        kallisto quant --bootstrap-samples 30 \\
                        --gtf ${gtf} \\
                        --threads 4 \\
						--rf-stranded \\
                        --index ${index} \\
						$endedness $unmapped \\
                        -o ${sample}
        """
        }
	}

/*
 * STEP 12b - Kallisto tximport to generate tpm and counts per gene
 */


    process kallisto_tximport {
      label 'low_memory'
      publishDir "${params.outdir}/kallisto", mode: 'copy'

      input:
      set val(name), file ('kallisto/*') from kallisto_tximport
      file gencode_csv from tx2gene_tximport.collect()

      output:
      file "${name}_kallisto_gene_tpm.csv" into kallisto_gene_tpm
      file "${name}_kallisto_gene_counts.csv" into kallisto_gene_counts
      file "${name}_kallisto_transcript_tpm.csv" into kallisto_transcript_tpm
      file "${name}_kallisto_transcript_counts.csv" into kallisto_transcript_counts

      script:
      """
      tximport.r NULL kallisto ${name} $gencode_csv
      """
    }

/*
 * STEP 12c - Kallisto merge for all samples
 */


    process kallisto_merge {
      label 'mid_memory'
      publishDir "${params.outdir}/kallisto", mode: 'copy'

      input:
      file gene_tpm_files from kallisto_gene_tpm.collect()
      file gene_count_files from kallisto_gene_counts.collect()
      file transcript_tpm_files from kallisto_transcript_tpm.collect()
      file transcript_count_files from kallisto_transcript_counts.collect()
      file gencode_csv from tx2gene_merge.collect()

      output:
      file "kallisto_merged*.csv" into kallisto_merged_ch
      file "*.rds"

      script:
      // First field is the gene/transcript ID
      gene_ids = "<(cut -f1 -d, ${gene_tpm_files[0]} | tail -n +2 | cat <(echo '${params.fc_group_features}') - )"
      transcript_ids = "<(cut -f1 -d, ${transcript_tpm_files[0]} | tail -n +2 | cat <(echo 'transcript_id') - )"

      // Second field is counts/TPM
      gene_tpm = gene_tpm_files.collect{f -> "<(cut -d, -f2 ${f})"}.join(" ")
      gene_counts = gene_count_files.collect{f -> "<(cut -d, -f2 ${f})"}.join(" ")
      transcript_tpm = transcript_tpm_files.collect{f -> "<(cut -d, -f2 ${f})"}.join(" ")
      transcript_counts = transcript_count_files.collect{f -> "<(cut -d, -f2 ${f})"}.join(" ")
      """
      paste -d, $gene_ids $gene_tpm > kallisto_merged_gene_tpm.csv
      paste -d, $gene_ids $gene_counts > kallisto_merged_gene_counts.csv
      paste -d, $transcript_ids $transcript_tpm > kallisto_merged_transcript_tpm.csv
      paste -d, $transcript_ids $transcript_counts > kallisto_merged_transcript_counts.csv
      se.r NULL kallisto_merged_gene_counts.csv kallisto_merged_gene_tpm.csv $gencode_csv
      se.r NULL kallisto_merged_transcript_counts.csv kallisto_merged_transcript_tpm.csv $gencode_csv
      """
    }


/*
 * STEP 13 - MultiQC
 */

/*
 * Creates a MultiQC Report
 * Following modules are incorporated:
 * Dupradar, biotype counts, edgeR MDS and similarity, preseq, resqc, featurecounts, 
 * rnaseqmetrics, STAR, cutadapt and fastqc on trimmed
 */

def nfcoreHeader() {
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\"${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\"${c_reset}
    ${c_purple}  hgx-rnaseq v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}


custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
  custom_runName = workflow.runName
}


// Header log info
//log.info nfcoreHeader() header when running pipeline, uncomment to trigger
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name'] = custom_runName ?: workflow.runName
summary['Reads'] = params.reads
summary['Data Type'] = params.singleEnd ? 'Single-End' : 'Paired-End'
if (params.genome) summary['Genome'] = params.genome

// Removed parameters

if (params.gtf) summary['GTF Annotation'] = params.gtf
if (params.gff) summary['GFF3 Annotation'] = params.gff
if (params.bed12) summary['BED Annotation'] = params.bed12
if (params.gencode) summary['GENCODE'] = params.gencode
summary['Remove Ribosomal RNA'] = params.removeRiboRNA
if (params.fc_group_features_type) summary['Biotype GTF field'] = biotype
summary['Save prefs'] = "Ref Genome: "+(params.saveReference ? 'Yes' : 'No')+" / Trimmed FastQ: "+(params.saveTrimmed ? 'Yes' : 'No')+" / Alignment intermediates: "+(params.saveAlignedIntermediates ? 'Yes' : 'No')
summary['Max Resources'] = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir'] = params.outdir
summary['Launch dir'] = workflow.launchDir
summary['Working dir'] = workflow.workDir
summary['Script dir'] = workflow.projectDir
summary['User'] = workflow.userName
if (workflow.profile == 'awsbatch') {
  summary['AWS Region']     = params.awsregion
  summary['AWS Queue']      = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
  summary['E-mail Address']    = params.email
  summary['E-mail on failure'] = params.email_on_fail
  summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"


 def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}


// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'hgx-rnaseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'hgx/rnaseq Workflow Summary'
    section_href: 'https://www.histogenex.com/'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
} 


process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    when:
    !params.skipMultiQC

    input:
    file multiqc_config from ch_multiqc_config
    file (fastqc:'fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('trimgalore/*') from trimgalore_results.collect().ifEmpty([])
    file ('alignment/*') from alignment_logs.collect().ifEmpty([])
    file ('rseqc/*') from rseqc_results.collect().ifEmpty([])
    file ('markDuplicates/*') from picard_results.collect().ifEmpty([])
    file ('RnaSeqMetrics/*') from rnaseqmetrics_results.collect().ifEmpty([])
    file ('preseq/*') from preseq_results.collect().ifEmpty([])
    file ('dupradar/*') from dupradar_results.collect().ifEmpty([])
    file ('featureCounts/*') from featureCounts_logs.collect().ifEmpty([])
    file ('featureCounts_biotype/*') from featureCounts_biotype.collect()
    file ('kallisto/*') from kallisto_logs.collect().ifEmpty([])
    file ('sortmerna/*') from sortmerna_logs.collect().ifEmpty([])
    file ('sample_correlation/*') from sample_correlation_results.collect().ifEmpty([])
    file workflow_summary from create_workflow_summary(summary)

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config \\
        -m custom_content -m picard -m preseq -m rseqc -m featureCounts -m star -m cutadapt -m sortmerna -m fastqc -m picard -m kallisto
    """
}


/*
 * Completion e-mail notification
 */


workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[hgx-rnaseq] Successful: $workflow.runName"
    if (skipped_poor_alignment.size() > 0) {
        subject = "[hgx-rnaseq] Partially Successful (${skipped_poor_alignment.size()} skipped): $workflow.runName"
    }
    if (!workflow.success) {
      subject = "[hgx-rnaseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if (workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['skipped_poor_alignment'] = skipped_poor_alignment
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success && !params.skipMultiQC) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[hgx-rnaseq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[hgx-rnaseq] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "/opt/conda/envs/hgx-rnaseq-1.0/", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
          if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[hgx-rnaseq] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, email_address ].execute() << email_txt
          log.info "[hgx-rnaseq] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = file("${output_d}/pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = file("${output_d}/pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (skipped_poor_alignment.size() > 0) {
        log.info "${c_purple}[hgx-rnaseq]${c_red} WARNING - ${skipped_poor_alignment.size()} samples skipped due to poor mapping percentages!${c_reset}"
    }
    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
        log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
        log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "${c_purple}[hgx-rnaseq]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[hgx-rnaseq]${c_red} Pipeline completed with errors${c_reset}"
    }

}
