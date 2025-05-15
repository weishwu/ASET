process aselux_index_genome {
    tag "aselux_index_genome"
    label "ASElux"

    input:
    tuple path(genome_fasta), path(genome_fai)
    path(gtf)

    output:
    tuple path("aselux_genome_genome.sa"), path("aselux_genome_gene.sa"), path("aselux_genome.annotation")

    script:
    """
    ASElux build --gtf ${gtf} --ref ${genome_fasta} --out aselux_genome 2>&1|tee >aselux_index_genome.log
    """
}

