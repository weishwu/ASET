process ase_hetsnp_add_genetype {
    tag "ase_hetsnp_add_genetype"
    label "pandas"

    input:
    path(ase_hetsnp)
    path(genes_bed)
    
    output:
    path("ASE_hetSNP_addExons_addGenetype.gz"), emit: ase_hetsnp

    script:
    """
    ase_hetsnp_add_genetype.py ${ase_hetsnp} ${genes_bed} ASE_hetSNP_addExons_addGenetype.gz 2>&1|tee >ASE_hetSNP_addGenetype.log
    """ 
}

