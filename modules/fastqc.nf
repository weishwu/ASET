process fastqc_raw {
    tag "fastqc_raw $sample"
    label "fastqc"

    publishDir "${params.dirs.results_dir}/qc/fastqc.raw_reads", mode: 'copy', pattern: "*.zip"
    publishDir "${params.dirs.results_dir}/qc/fastqc.raw_reads/", mode: 'copy', pattern: "*.html"

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path("*_fastqc.zip"), emit: fastqc_zip
    path("*_fastqc.html"), emit: fastqc_html

    script:
    """
    fastqc -o . ${read1}
    fastqc -o . ${read2}
    """
}

process fastqc_trimmed {
    tag "fastqc_trimmed $sample"
    label "fastqc"

    publishDir "${params.dirs.results_dir}/qc/fastqc.trimmed_reads", mode: 'copy', pattern: "*.zip"
    publishDir "${params.dirs.results_dir}/qc/fastqc.trimmed_reads/", mode: 'copy', pattern: "*.html"

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path("*_fastqc.zip"), emit: fastqc_zip
    path("*_fastqc.html"), emit: fastqc_html

    script:
    """
    fastqc -o . ${read1}
    fastqc -o . ${read2}
    """
}

