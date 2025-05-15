process star_genome_index {
    tag "star_genome_index"
    label "STAR"

    input: 
    tuple path(genome_fasta), path(genome_fai)
    path(gtf)

    output:
    path("STAR_genome"), emit: star_genome

    script:
    """
    overhang=\$(echo ${params.data.rna_readlen} | awk '{print \$1-1}')
    STAR --runThreadN ${task.cpus} \
      --runMode genomeGenerate \
      --genomeDir STAR_genome \
      --genomeFastaFiles ${genome_fasta} \
      --sjdbGTFfile ${gtf} \
      --sjdbOverhang \${overhang} \
      --limitGenomeGenerateRAM ${task.memory.toGiga()}000000000 \
    2>&1|tee >STAR_genomeGenerate.log 2>&1
    """ 
}



