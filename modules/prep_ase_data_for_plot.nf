process prep_ase_data_for_plot {
    tag "prep_ase_data_for_plot"
    label "rplot"

    publishDir "${params.dirs.results_dir}/ase_data/", mode: 'copy', pattern: "*Rds"

    input:
    path(ase_hetsnp)
    path(merged_exon_bed)
    
    output:
    path("ase_data_for_plot.Rds"), emit: ase_data_for_plot

    script:
    """
    prep_ase_data_for_plot.R ${ase_hetsnp} ${merged_exon_bed} ${params.data.ase_depth_min} ase_data_for_plot.Rds \
    2>&1|tee >prep_ase_data_for_plot.log
    """ 
}

