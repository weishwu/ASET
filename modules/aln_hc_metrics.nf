process aln_hc_metrics {
    tag "aln_hc_metrics: ${sample}"
    label "gatk"

    publishDir "${params.dirs.results_dir}/qc/aln_hc_metrics/", mode: 'copy', pattern: "*txt"

    input:
    tuple val(sample), path(aln_hc_bam)
    path(gtf_refflat)
    path(rRNA_tRNA_list)
    
    output:
    path("${sample}.hc.bam.RNAseqMetrics.txt")

    script:
    """
    #export JAVA_OPTIONS=-Xmx${task.memory.toGiga()}g
    gatk CollectRnaSeqMetrics \
        --java-options "-Xmx${task.memory.toGiga()}g" \
        --REF_FLAT ${gtf_refflat} \
        --RIBOSOMAL_INTERVALS ${rRNA_tRNA_list} \
        --STRAND_SPECIFICITY ${params.data.rna_strandedness} \
        -I ${aln_hc_bam} \
        -O ${sample}.hc.bam.RNAseqMetrics.txt \
        --TMP_DIR . \
    2>&1|tee > ${sample}.hc.bam.RNAseqMetrics.log
    """ 
}
