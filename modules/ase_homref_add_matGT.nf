process ase_homref_add_matGT {
    tag "ase_homref_add_matGT"
    label "pandas"
    label "py_parallel"

    input:
    path(ase_count_homref)
    path(mat_vcf_sheet)
    path(mat_vcf_files)
    
    output:
    path("ASE_homRef_MatGT.gz"), emit: ase_homref_with_matGT

    script:
    """
    ase_homref_add_matGT.py ${ase_count_homref} ${mat_vcf_sheet} ASE_homRef_MatGT.gz ${task.cpus} 2>&1|tee >ase_homref_add_matGT.log
    """ 
}

