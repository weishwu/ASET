process gmap_index_genome {
    tag "gmap_index_genome"
    label "GSNAP"

    input: 
    tuple path(genome_fasta), path(genome_fai)

    output:
    path("gmap_genome"), emit: gmap_genome

    script:
    """
    # -s numeric-alpha
    gmap_build \
      -t ${task.cpus} \
      -D gmap_genome \
      -d gmap_genome \
      ${genome_fasta} \
    >gmap_index_genome.log 2>&1
    """ 
}

