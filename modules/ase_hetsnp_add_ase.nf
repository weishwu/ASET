process ase_hetsnp_add_ase {
    tag "ase_hetsnp_add_ase"
    label "pandas"

    input:
    path(ase_hetsnp)
    
    output:
    path("ASE_hetSNP_addExons_addmatcontam_addphasing_addase.gz"), emit: ase_hetsnp

    script:
    """
    ase_hetsnp_add_ase.py ${ase_hetsnp} ASE_hetSNP_addExons_addmatcontam_addphasing_addase.gz 2>&1|tee >ASE_hetSNP_addExons_addmatcontam_addphasing_addase.log
    """ 
}

