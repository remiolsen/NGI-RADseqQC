#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                                     RADQC Pipeline 
========================================================================================
 QC pipeline for NGI genotyping-by-sequencing (RAD-seq) data.
 @Authors
 Remi-Andre Olsen <remi-andre.olsen@scilifelab.se> 
----------------------------------------------------------------------------------------
 Basic command:
 $ nextflow main.nf
 
 Pipeline parameters
 --fastqpath
 --trim-adapters
 --trim-truncate
 --enz
 --outdir
 --project
 --clusterOptions
 --small-m
 --big-m
 --small-n

----------------------------------------------------------------------------------------
*/
/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 0.2

params.fastqpath = "${baseDir}/example/P*/*/*R{1,2}_001.fastq.gz"
params.trim_adapters = "${baseDir}/resources/nextera_linkers.txt"
params.trim_truncate = 100
params.trimmomatic_home = null
params.enz = "ecoRI"
params.outdir = "$PWD"
params.project = "b2013064"
params.small_m = 3
params.big_m = 2
params.small_n = 1
params.clusterOptions = null
if (params.trimmomatic_home) {
    env = {TRIMMOMATIC_HOME = params.trimmomatic_home}
}

log.info "###"
log.info "###    NGI RADQC pipeline v${version}"
log.info "###"
params.each { key, value ->
    log.info "$key = $value"
}
log.info "denovo_map.pl to use parameters: -m ${params.small_m} -M ${params.big_m} -n ${params.small_n}"

read_files_trim = Channel
    // TODO: replace with a better pattern
    .fromFilePairs( "${params.fastqpath}", size: -1 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.fastqpath}" }




process trimmomatic {
    tag "$name"
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'
    
    input:
    set val(name), file(reads) from read_files_trim

    output:
    //consider to output unpaired as well to jellyfish
    set val(name), file("*P.fastq.gz") into read_files_flash, read_files_jellyfish, readP_files_concat
    set val(name), file("*U.fastq.gz") into readU_files_concat
    file "*_trim.out"

    script:
    """
    java -jar \$TRIMMOMATIC_HOME/trimmomatic.jar PE \
    -threads 1 \
    -trimlog ${name}_trim.log \
    -baseout ${name}.fastq.gz \
    -phred33 $reads \
    ILLUMINACLIP:${params.trim_adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 \
    MINLEN:${params.trim_truncate} CROP:${params.trim_truncate} 2> ${name}_trim.out 
    """

}


process flash {
    tag "$name"
    publishDir "${params.outdir}/flash", mode: 'copy'
    
    input:
    set val(name), file(reads) from read_files_flash
    
    output:
    file "*_flash.txt" into flash_results

    script:
    """
    flash -t ${task.cpus} -M ${params.trim_truncate} -c $reads 2> ${name}_flash.txt > /dev/null
    """

}

process concat_reads {
    tag "$name"
    publishDir "${params.outdir}/concatenated_reads", mode: 'copy'

    input:
    set val(name), file(readsP) from readP_files_concat
    set val(name), file(readsU) from readU_files_concat

    output:
    set val(name), file("*.merged.fastq.gz") into merged_reads 

    script:
    """
    cat $readsP $readsU > ${name}.merged.fastq.gz
    """

} 


process process_radtags {
    tag "$name"
    publishDir "${params.outdir}/process_radtags", mode: 'copy',
    saveAs: {it == 'process_radtags.log' ? name+'_process_radtags.log' : it}

    input:
    set val(name), file(reads) from merged_reads

    output:
    file "*process_radtags.log"
    file "*.merged.fq.gz" into processed_reads, read_files_fastqc

    script:
    """
    process_radtags -i gzfastq -f $reads -e ${params.enz} -c -q -r -o .
    """
}

/*
    If you are here, you should already have run FastQC on the raw input
*/
process fastqc {
    tag "$reads"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}"

    script:
    """
    fastqc -q $reads
    """
}


process denovo_stacks {
    publishDir "${params.outdir}/denovo_stacks", mode: 'copy'
    
    input:
    file(reads) from processed_reads.toList()

    output:
    file "*.tsv.gz"
    file "*.tsv"
    file "denovo_map.log" into denovo_log
    file "*.populations.log"

    script:
    s_string = ""
    reads.sort(); reads.each {s_string = s_string + "-s $it "}
    """
    denovo_map.pl $s_string -o . -m ${params.small_m} -M ${params.big_m} -n ${params.small_n} -S -b 1 -n 2 -T ${task.cpus}
    """
}
