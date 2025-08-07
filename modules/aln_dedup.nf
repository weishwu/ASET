process aln_dedup {
    tag "aln dedup: ${sample}"
    label "gatk"

    publishDir "${params.dirs.results_dir}/aln/03.aln_hc_dedup/", mode: 'copy', pattern: "*ba?"
    publishDir "${params.dirs.results_dir}/aln/03.aln_hc_dedup/", mode: 'copy', pattern: "*metrics.txt"
    publishDir "${params.dirs.results_dir}/aln/03.aln_hc_dedup/", mode: 'copy', pattern: "*flagstat"

    input:
    tuple val(sample), path(aln_hc_bam)
    
    output:
    tuple val(sample), path("${sample}.aln.hc.dedup.bam"), emit: aln_hc_dedup_bam
    path("aln.hc_filt.dedup.${sample}.metrics.txt"), emit: aln_hc_dedup_metrics
    tuple val(sample), path("${sample}.aln.hc.dedup.bai")

    script:
    """
    #export JAVA_OPTIONS="-Xmx${task.memory.toGiga()}g"
    gatk MarkDuplicates \
        --java-options "-Xmx${task.memory.toGiga()}g" \
        -I ${aln_hc_bam} \
        -O ${sample}.aln.hc.dedup.bam \
        --METRICS_FILE aln.hc_filt.dedup.${sample}.metrics.txt \
        --VALIDATION_STRINGENCY SILENT \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --ASSUME_SORT_ORDER coordinate \
        --CREATE_INDEX true \
        --REMOVE_DUPLICATES true \
        --TMP_DIR ./ \
    2>&1 |tee > aln.hc_filt.dedup.${sample}.log
    """ 
}

