process prep_diploid_genome {
    tag "prep_diploid_genome: ${sample}"
    label "vcf2diploid"

    input: 
    path(genome)
    tuple val(sample), path(snps)

    output:
    tuple val(sample), path("${sample}_HAP1/*.fasta"), path("${sample}_HAP2/*.fasta"), emit: diploid_fa 
    tuple val(sample), path("${sample}_HAP1/*.chain"), path("${sample}_HAP2/*.chain"), emit: diploid_chain

    script:
    """
    java -Xmx${task.memory.toGiga()}g -jar vcf2diploid.jar -id ${sample} -chr ${genome} -vcf ${snps} > ${sample}.vcf2diploid.logfile.txt
    """ 
}


