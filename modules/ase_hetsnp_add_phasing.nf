process ase_hetsnp_add_phasing {
    tag "ase_hetsnp_add_phasing"
    label "pandas"

    input:
    path(ase_hetsnp)
    path(phase_info_sheet)
    path(phase_info_files)
    
    output:
    path("ASE_on_hetSNP.gz"), emit: ase_hetsnp

    script:
    """
    ase_hetsnp_add_phasing.py ${ase_hetsnp} ${phase_info_sheet} ASE_on_hetSNP.gz 2>&1|tee >ASE_hetSNP_addExons_addmatcontam_addphasing.log
    """ 
}

