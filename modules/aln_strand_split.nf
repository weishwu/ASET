process aln_strand_split {
    tag "aln strand split: ${sample}"
    label "ngsutils"

    publishDir "${params.dirs.results_dir}/aln/04.aln_hc_dedup_strandsplit/", mode: 'copy', pattern: "*ba?"
    publishDir "${params.dirs.results_dir}/aln/03.aln_hc_dedup/", mode: 'copy', pattern: "*flagstat"

    input:
    tuple val(sample), path(dedup_bam)
    
    output:
    tuple val(sample), path("${sample}.aln.hc.dedup.sense.bam"), path("${sample}.aln.hc.dedup.antisense.bam"), path("${sample}.aln.hc.dedup.sense.bam.bai"), path("${sample}.aln.hc.dedup.antisense.bam.bai"), emit: aln_strand_split_bam
    path("${sample}.aln.hc.dedup.bam.flagstat"), emit: aln_hc_dedup_flagstat
 
    script:
    """
    half_cpus=`echo ${task.cpus} | awk '{print \$1/2}'`
    samtools flagstat -@${task.cpus} ${dedup_bam} > ${sample}.aln.hc.dedup.bam.flagstat
    samtools view -@\${half_cpus} -h ${dedup_bam} | \
       awk '{if ((\$1 ~ "^@") || (\$2 == 83) || (\$2 == 163)) print \$0}' | \
       samtools view -@\${half_cpus} -b >${sample}.aln.hc.dedup.r2_sense.bam
    samtools view -@\${half_cpus} -h ${dedup_bam} | \
       awk '{if ((\$1 ~ "^@") || (\$2 == 99) || (\$2 == 147)) print \$0}' | \
       samtools view -@\${half_cpus} -b >${sample}.aln.hc.dedup.r2_antisense.bam
    if [[ ${params.data.rna_strandedness} == 'SECOND_READ_TRANSCRIPTION_STRAND' ]]; then
       mv ${sample}.aln.hc.dedup.r2_sense.bam ${sample}.aln.hc.dedup.sense.bam
       mv ${sample}.aln.hc.dedup.r2_antisense.bam ${sample}.aln.hc.dedup.antisense.bam
    elif [[ ${params.data.rna_strandedness} == 'FIRST_READ_TRANSCRIPTION_STRAND' ]]; then
       mv ${sample}.aln.hc.dedup.r2_antisense.bam ${sample}.aln.hc.dedup.sense.bam
       mv ${sample}.aln.hc.dedup.r2_sense.bam ${sample}.aln.hc.dedup.antisense.bam
    else
       echo "rna_strandedness parameter has to be either 'SECOND_READ_TRANSCRIPTION_STRAND' or 'FIRST_READ_TRANSCRIPTION_STRAND'!" && exit
    fi
    samtools index ${sample}.aln.hc.dedup.sense.bam
    samtools index ${sample}.aln.hc.dedup.antisense.bam
    """ 
}

