process ase_count_split {
    tag "ase_count_split: ${sample}"
    label "pandas"

    input:
    tuple val(sample), path(sense_ase), path(antisense_ase), path(snps)

    output:
    path("*.aln.hc.dedup.ASE.hetSNP.txt"), emit: ase_count_hetsnp_txt
    path("*.aln.hc.dedup.ASE.homSNP.txt"), emit: ase_count_homsnp_txt
    path("*.aln.hc.dedup.ASE.homRef.txt"), emit: ase_count_homref_txt

    script:
    """
    ase_count_split.py ${sample} ${snps} ${sense_ase} ${antisense_ase} ${sample}.aln.hc.dedup.ASE.hetSNP.txt ${sample}.aln.hc.dedup.ASE.homSNP.txt ${sample}.aln.hc.dedup.ASE.homRef.txt \
    2>&1|tee > ase_count_split.${sample}.log
    """ 
}


