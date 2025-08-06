process gtf_to_df {
    tag "gtf_to_df"
    label "pandas"
    label "py_parallel"

    publishDir "${params.dirs.results_dir}/ref_data/", mode: 'copy', pattern: "*.txt"

    input:
    path(gtf)

    output:
    path("gtf_flat.txt"), emit: gtf_to_df_txt

    script:
    """
    flat_gtf.py ${gtf} gtf_flat.txt ${task.cpus} \
    2>&1|tee > flat_gtf.log
    """ 
}

