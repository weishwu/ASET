process aln_ase_count {
    tag "aln ase_count: ${sample}"
    label "gatk"

    input:
    tuple path(genome_fa), path(genome_fai)
    path(genome_dict)
    tuple val(sample), path(sense_bam), path(antisense_bam), path(sense_bai), path(antisense_bai), path(snps)

    output:
    tuple val(sample), path("${sample}.aln.hc.dedup.sense.ASE.txt"), path("${sample}.aln.hc.dedup.antisense.ASE.txt")

    script:
    """
    #export JAVA_OPTIONS=-Xmx${task.memory.toGiga()}g
    
    gunzip -c ${snps} | awk '{OFS="\\t"; if (\$1 ~ "^#") {print \$0} else {\$9="GT";\$10="0/1"; print \$0}}' > ${sample}.pseudoHet.vcf
    gatk SortVcf -I ${sample}.pseudoHet.vcf -O ${sample}.pseudoHet.vcf.gz --CREATE_INDEX true

    for strand in sense antisense
    do 
    gatk ASEReadCounter ${params.tool_parameters.asereadcounter_flags} \
        --java-options "-Xmx${task.memory.toGiga()}g" \
        -R ${genome_fa} \
        -I ${sample}.aln.hc.dedup.\${strand}.bam \
        -V ${sample}.pseudoHet.vcf.gz \
        -O ${sample}.aln.hc.dedup.\${strand}.ASE.txt \
        --output-format TABLE \
        --tmp-dir . \
    2>&1|tee > ase_read_counter.${sample}.\${strand}.log
    done
    rm -f ${sample}.pseudoHet.vcf ${sample}.pseudoHet.vcf.gz
    """ 
}
