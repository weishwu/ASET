process genome_link {
    tag "genome_link"
    label "ngsutils"

    input:
    path(genome_fa)

    output:
    tuple path("genome.fa"), path("genome.fa.fai")

    script:
    """
    ln -s ${genome_fa} genome.fa
    samtools faidx genome.fa
    """
}

process genome_dict {
    tag "genome_dict"
    label "gatk"

    input:
    tuple path(genome_fa), path(genome_fai)

    output:
    path("genome.dict")

    script:
    """
    gatk CreateSequenceDictionary \
        -R ${genome_fa} \
        -O genome.dict
    """
}

