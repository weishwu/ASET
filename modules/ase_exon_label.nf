process ase_exon_label {
    tag "ase_exon_label"
    label "pandas"

    input:
    path(ase_count_hetsnp)
    path(ase_x_exons)
    path(ase_x_genes)
    
    output:
    path("ASE_hetSNP_addExons.gz"), emit: ase_exon_label

    script:
    """
    ase_exon_label.py ${ase_count_hetsnp} ${ase_x_exons} ${ase_x_genes} ASE_hetSNP_addExons.gz
    """ 
}

