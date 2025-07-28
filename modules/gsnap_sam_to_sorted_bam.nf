process gsnap_sam_to_sorted_bam {
    tag "gsnap_sam_to_sorted_bam: ${sample}"
    label "ngsutils"

    publishDir "${params.dirs.results_dir}/aln/01.raw_aln", mode: 'copy', pattern: "*ba?"

    input:
    tuple val(sample), path(sam)

    output:
    tuple val(sample), path("${sample}.gsnap.sorted.bam")

    script:
    """
    half_cpus=`echo ${task.cpus} | awk '{print \$1/2}'`
    samtools view -@\${half_cpus} -bS ${sam} | samtools sort -@\${half_cpus} -O BAM -o ${sample}.gsnap.sorted.bam
    samtools index ${sample}.gsnap.sorted.bam
    """
}
