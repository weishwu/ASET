process extract_exons_genes {
    tag "extract_exons_genes"
    label "pandas"

    publishDir "${params.dirs.results_dir}/ref_data/", mode: 'copy', pattern: "*bed"

    input:
    path(gtf_flat)
    
    output:
    path("ref_gtf_exons.bed"), emit: exons_bed
    path("ref_gtf_genes.bed"), emit: genes_bed

    script:
    """
    extract_exons_genes.py ${gtf_flat} ref_gtf_exons.bed ref_gtf_genes.bed \
    2>&1|tee > extract_exons_genes.log
    """ 
}


