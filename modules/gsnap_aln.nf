process gsnap_aln {
    tag "gsnap_aln: ${sample}"
    label "GSNAP"

    input: 
    path(gtf_index)
    tuple val(sample), path(trimmed_read1), path(trimmed_read2), path(genome_index), path(snps)

    output:
    tuple val(sample), path("${sample}.gsnap.sam"), emit: aln_sam
    path("gsnap_aln.${sample}.log"), emit: aln_log

    script:
    """
    gsnap \
       --read-group-id ${sample} \
       --read-group-name ${sample} \
       -D ${genome_index} \
       -d gmap_genome \
       -s ${gtf_index} \
       -v ${sample} \
       -t ${task.cpus} \
       --gunzip \
       --format sam \
       -o ${sample}.gsnap.sam \
       ${params.tool_parameters.gsnap} \
       ${trimmed_read1} \
       ${trimmed_read2} \
    >gsnap_aln.${sample}.log 2>&1
    """ 
}


