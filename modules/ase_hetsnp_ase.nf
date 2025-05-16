process ase_hetsnp_ase {
    tag "ase_hetsnp_ase"
    label "pandas"

    input:
    path(ase_hetsnp)
    path(ase_homsnp)
    path(ase_homref)
    
    output:
    path("ASE_hetSNP_addExons_addcontam.gz"), emit: ase_hetsnp

    script:
    """
    ase_hetsnp_ase.py ${ase_hetsnp} ${ase_homsnp} ${ase_homref} ${params.data.ase_depth_min} ASE_hetSNP_addExons_addcontam.gz  2>&1|tee >ASE_hetSNP_addExons_addcontam.log
    """ 
}

