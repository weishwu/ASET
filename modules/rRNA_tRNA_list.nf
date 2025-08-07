process rRNA_tRNA_list {
    tag "rRNA_tRNA_list"
    label "gatk"

    input:
    path(bed)
    path(genome_dict)

    output:
    path("rRNA_tRNA_list.interval_list")

    script:
    """
    #export JAVA_OPTIONS=-Xmx${task.memory.toGiga()}g
    gatk BedToIntervalList \
        --java-options "-Xmx${task.memory.toGiga()}g" \
        -I ${bed} \
        -SD ${genome_dict} \
        -O rRNA_tRNA_list.interval_list
    """
}

