process aselux_ase {
    tag "aselux: ${sample}"
    label "ASElux"

    input:
    tuple path(aselux_genome), path(aselux_gene), path(aselux_annotation)
    tuple val(sample), path(trimmed_read1), path(trimmed_read2), path(snps)

    output:
    tuple val(sample), path("aselux.${sample}.txt")

    script:
    """
    gunzip -c ${snps} > ${sample}.snps.vcf
    gunzip -c ${trimmed_read1} > ${sample}.R1.fastq
    gunzip -c ${trimmed_read2} > ${sample}.R2.fastq
    ASElux align --nthread 1 --fq --pe --readLen ${params.data.rna_readlen} --index aselux_genome --vcf ${sample}.snps.vcf --seqFiles ${sample}.R1.fastq ${sample}.R2.fastq --out aselux.${sample}.txt 2>&1|tee >aselux.${sample}.log
    rm ${sample}.snps.vcf ${sample}.R1.fastq ${sample}.R2.fastq
    """
}

