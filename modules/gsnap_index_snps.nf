process gsnap_index_snps {
    tag "gsnap_index_snps: ${sample}"
    label "GSNAP"

    input: 
    path(gmap_genome)
    tuple val(sample), path(snps)

    output:
    tuple val(sample), path("gmap_genome"), path("gmap_genome/gmap_genome/gmap_genome.ref153offsets64meta.${sample}"), emit: gsnap_snps

    script:
    """
gunzip -c ${snps} | vcf_iit > ${sample}_snps.txt
cat ${sample}_snps.txt | iit_store -o ${sample}_snps
snpindex -D ${gmap_genome} -d gmap_genome -v ${sample} ${sample}_snps.iit
    """ 
}


