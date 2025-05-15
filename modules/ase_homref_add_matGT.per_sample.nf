process ase_homref_add_matGT {
    tag "ase_homref_add_matGT: ${sample}"
    label "pandas"

    input:
    tuple val(sample), path(ase_count_homref), path(mat_snp)
    
    output:
    tuple val(sample), path("ASE_homRef_MatGT.${sample}.txt"), emit: ase_homref_with_matGT

    script:
    """
    ase_homref_add_matGT.py ${ase_count_homref} ${mat_snp} ASE_homRef_MatGT.${sample}.txt
    """ 
}

