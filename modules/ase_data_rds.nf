process ase_data_rds {
    tag "ase_data_rds"
    label "rplot"

    publishDir "${params.dirs.results_dir}/ase_data/", mode: 'copy', pattern: "*Rds"
    publishDir "${params.dirs.results_dir}/ase_data/", mode: 'copy', pattern: "ASE_on_hetSNPs.txt.gz"

    input:
    path(ase_hetsnp)
    path(merged_exon_bed)
    
    output:
    path("ase_data.Rds"), emit: ase_data_rds
    path("ASE_on_hetSNPs.txt.gz"), emit: ase_data_txt

    script:
    """
    ase_data_rds.R ${ase_hetsnp} ${merged_exon_bed} ase_data.Rds
    2>&1|tee >ase_data_rds.log
    cp ${ase_hetsnp} ASE_on_hetSNPs.txt.gz
    """ 
}

