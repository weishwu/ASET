process star_nmaskedgenome_aln {
    tag "star_nmaskedgenome_aln: ${sample}"
    label "STAR"

    publishDir "${params.dirs.results_dir}/aln/01.raw_aln/", mode: 'copy', pattern: "*ba?"
    publishDir "${params.dirs.results_dir}/aln/01.raw_aln/", mode: 'copy', pattern: "*.Log.final.out"

    input:
    tuple val(sample), path(star_genome), path(trimmed_read1), path(trimmed_read2)

    output:
    tuple val(sample), path("${sample}.Aligned.sortedByCoord.out.bam"), emit: aln_bam
    path("${sample}.Log.final.out"), emit: aln_log

    script:
    """
    STAR \
        --genomeDir ${star_genome} \
        --runThreadN ${task.cpus} \
        --readFilesIn ${trimmed_read1} ${trimmed_read2} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${sample}. \
        --outSAMattrRGline ID:${sample} SM:${sample} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI NM MD AS nM jM jI XS \
        ${params.tool_parameters.star} \
    2>&1|tee >star_aln.${sample}.log 
    """ 
}

