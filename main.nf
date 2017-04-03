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
 --fastqdir
 --trim
 --trim-adapters
 --trim-truncate
 --enz
 --outdir
 
----------------------------------------------------------------------------------------
*/
/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 0.1

params.fastqdir = "${baseDir}/example/"
params.trim_adapters = "${baseDir}/resources/nextera_linkers.txt"
params.trim = true
params.trim_truncate = 100
params.enz = "ecoRI"
params.outdir = "$PWD"



log.info "### RADQC pipeline v${version}"
log.info "fastqdir = ${params.fastqdir}"
log.info "trim = ${params.trim}"
log.info "trim_adapters = ${params.trim_adapters}"
log.info "trim_truncate = ${params.trim_truncate}"
log.info "enz = ${params.enz}"
log.info "outdir = ${params.outdir}"


/*
 * Always start with fastqc 
 */
Channel
    // TODO: replace with a better pattern
    .fromFilePairs( "${params.fastqdir}/P*/*/*R{1,2}_001.fastq.gz", size: -1 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.fastqdir}" }
    .into { read_files_fastqc; read_files_trim }


process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}"

    script:
    """
    fastqc -q $reads
    """
}

if(params.trim) {

    process trimmomatic {
        tag "$name"
        publishDir  "${params.outdir}/trimmed_reads", mode: 'copy'
        
        input:
        set val(name), file(reads) from read_files_trim

        output:
        //output unpaired as well
        set val(name), file("*P.fastq.gz") into read_files_flash, read_files_jellyfish, read_files_concat
        file "*_trim.out"

        when:
        params.trim == true

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

}
else {
    read_files = Channel.fromFilePairs( "${params.fastqdir}/P*/*/*R{1,2}_001.fastq.gz", size: -1 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.fastqdir}" }
        .into { read_files_flash; read_files_jellyfish; read_files_concat }

}


process flash {
    tag "$name"
    publishDir "${params.outdir}/flash", mode: 'copy'
    
    input:
    set val(name), file(reads) from read_files_flash
    
    output:
    file "*_flash.txt" into flash_results

    // TODO: specify threads
    script:
    """
    flash -t ${task.cpus} -M ${params.trim_truncate} -c $reads 2> ${name}_flash.txt > /dev/null
    """

}

process jellyfish {
    tag "$name"
    publishDir "${params.outdir}/jellyfish", mode: 'copy'
    
    input:
    set val(name), file(reads) from read_files_jellyfish

    output:
    file "*.hist" into jellyfish_results

    //TODO: specify threads & mem
    script:
    """
    cat $reads | gzip -c -d | jellyfish count -o ${name}.jf -m 25 -s 1000M -t ${task.cpus} -C /dev/fd/0
    jellyfish histo -o ${name}.hist -f ${name}.jf
    """

}

process concat_reads {
    tag "$name"
    publishDir  "${params.outdir}/concatenated_reads", mode: 'copy'

    input:
    set val(name), file(reads) from read_files_concat

    output:
    set val(name), file("*.merged.fastq.gz") into merged_reads 

    script:
    """
    cat $reads > ${name}.merged.fastq.gz
    """

} 


process process_radtags {
    tag "$name"
    publishDir "${params.outdir}/process_radtags", mode: 'copy'

    input:
    set val(name), file(reads) from merged_reads

    output:
    stdout into process_logs 
    file "*.merged.fq.gz" into processed_reads

    script:
    """
    process_radtags -i gzfastq -f $reads -e ${params.enz} -c -q -r -o .
    cat process_radtags.log
    """
}
process_logs
  .collectFile(name: file("process_radtags.log"), storeDir: "process_radtags")
  .println { "Result saved to file: $it" }


process denovo_stacks {
    publishDir "${params.outdir}/denovo_stacks", mode: 'copy'
    
    input:
    file(reads: 'process_radtags/*.fq.gz') from processed_reads.toList()

    output:
    file "*.tsv.gz"
    file "*.tsv"
    file "denovo_map.log" into denovo_log
    file "*.populations.log"

    script:
    s_string = ""
    reads.each {s_string = s_string + "-s $it "}
    """
    denovo_map.pl $s_string -o . -m 3 -M 2 -N 4 -S -b 1 -n 2 -T ${task.cpus}
    """
}
