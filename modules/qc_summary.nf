process qc_summary_from_fastq {
    tag "qc_summary_from_fastq"
    label "pandas"
    label "py_parallel"

    publishDir "${params.dirs.results_dir}/qc/", mode: 'copy', pattern: "*txt"

    when:
    (params.data.routine == 'from_fastq') && (params.tool_parameters.mapper != 'ASElux')

    input:
    path(ase_hetsnp)
    path(sample_sheet)
    path(trimmomatic_logs)
    path(aln_logs)
    path(hc_flagstat)
    path(dedup_flagstat)
    path(hc_metrics)
    path(snps)

    output:
    path("QC_summary.txt"), emit: ase_hetsnp

    script:
    """
    if [[ (${params.tool_parameters.mapper} == 'STAR_WASP') || (${params.tool_parameters.mapper} == 'STAR_NMASK') ]]; then
      qc_summary_from_fastq_star_aln.py ${ase_hetsnp} ${sample_sheet} QC_summary.txt ${params.data.ase_depth_min} ${task.cpus} 2>&1|tee >qc_summary.log
    elif [[ ${params.tool_parameters.mapper} == 'GSNAP' ]]; then
      qc_summary_from_fastq_gsnap_aln.py ${ase_hetsnp} ${sample_sheet} QC_summary.txt ${params.data.ase_depth_min} ${task.cpus} 2>&1|tee >qc_summary.log
    fi
    """ 
}

process qc_summary_from_bam {
    tag "qc_summary_from_bam"
    label "pandas"
    label "py_parallel"

    publishDir "${params.dirs.results_dir}/qc/", mode: 'copy', pattern: "*txt"

    when:
    params.data.routine == 'from_bam'

    input:
    path(ase_hetsnp)
    path(sample_sheet)
    path(hc_flagstat)
    path(dedup_flagstat)
    path(hc_metrics)
    path(snps)

    output:
    path("QC_summary.txt"), emit: ase_hetsnp

    script:
    """
    qc_summary_from_bam.py ${ase_hetsnp} ${sample_sheet} QC_summary.txt ${params.data.ase_depth_min} ${task.cpus} \
    2>&1|tee >qc_summary.log
    """ 
}

process qc_summary_from_fastq_aselux {
    tag "qc_summary_from_fastq_aselux"
    label "pandas"
    label "py_parallel"

    publishDir "${params.dirs.results_dir}/qc/", mode: 'copy', pattern: "*txt"

    when:
    (params.data.routine == 'from_fastq') && (params.tool_parameters.mapper == 'ASElux')

    input:
    path(ase_hetsnp)
    path(sample_sheet)
    path(trimmomatic_logs)
    path(snps)

    output:
    path("QC_summary.txt"), emit: ase_hetsnp

    script:
    """
    qc_summary_from_fastq_aselux.py ${ase_hetsnp} ${sample_sheet} QC_summary.txt ${params.data.ase_depth_min} ${task.cpus} 2>&1|tee >qc_summary.log
    """ 
}
