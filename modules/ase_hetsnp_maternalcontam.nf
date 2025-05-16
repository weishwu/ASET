process ase_hetsnp_maternalcontam {
    tag "ase_hetsnp_maternalcontam"
    label "pandas"

    input:
    path(ase_hetsnp)
    path(ase_homref_exonic)
    
    output:
    path("ASE_hetSNP_addExons_addmatcontam.gz"), emit: ase_hetsnp

    script:
    """
    ase_hetsnp_maternalcontam.py ${ase_hetsnp} ${ase_homref_exonic} ${params.data.ase_depth_min} ASE_hetSNP_addExons_addmatcontam.gz 2>&1|tee >ASE_hetSNP_addExons_addmatcontam.log
    """ 
}

