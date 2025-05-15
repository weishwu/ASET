process star_nmaskedgenome_index {
    tag "star_nmaskedgenome_index: ${sample}"
    label "STAR"

    input: 
    path(gtf)
    tuple val(sample), path(genome_fasta)
  
    output:
    tuple val(sample), path("STAR_genome.${sample}")

    script:
    """
    overhang=\$(echo ${params.data.rna_readlen} | awk '{print \$1-1}')
    STAR --runThreadN ${task.cpus} \
      --runMode genomeGenerate \
      --genomeDir STAR_genome.${sample} \
      --genomeFastaFiles ${genome_fasta} \
      --sjdbGTFfile ${gtf} \
      --sjdbOverhang \${overhang} \
      --limitGenomeGenerateRAM ${task.memory.toGiga()}000000000 \
    2>&1|tee >star_genome_generate.${sample}.log
    """
}

