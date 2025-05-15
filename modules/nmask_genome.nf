process nmask_genome {
    tag "nmask_genome: ${sample}"
    label "bedtools"

    input: 
    tuple path(genome_fasta), path(genome_fai)
    tuple val(sample), path(snps)

    output:
    tuple val(sample), path("genome_Nmasked.${sample}.fa")

    script:
    """
    bedtools maskfasta \
       -fi ${genome_fasta} \
       -bed ${snps} \
       -fo genome_Nmasked.${sample}.fa 
    """ 
}

