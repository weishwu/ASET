process multiqc_from_fastq {
  tag "multiqc_from_fastq"
  label "multiqc"

  publishDir "${params.dirs.results_dir}/qc/", mode: 'copy', pattern: "*data"
  publishDir "${params.dirs.results_dir}/qc/", mode: 'copy', pattern: "*.html"

  input:
  path(fastqc_raw_zip)
  path(fastqc_raw_html)
  path(trimmomatic_log)
  path(fastqc_trimmed_zip)
  path(fastqc_trimmed_html)
  path(star_aln_log)
  path(aln_hc_metrics)
  path(aln_hc_dedup_metrics)

  output:
  path("*data")
  path("*html")

  script:
  """
  mkdir pre_trim
  mv ${fastqc_raw_zip} pre_trim
  mv ${fastqc_raw_html} pre_trim
  multiqc --filename pre_trim.multiqc ./pre_trim 2>&1|tee > multiqc.pre_trim.log
  multiqc --filename post_trim.multiqc --ignore pre_trim . 2>&1|tee > multiqc.post_trim.log
  """
}

process multiqc_from_bam {
  tag "multiqc_from_bam"
  label "multiqc"

  publishDir "${params.dirs.results_dir}/qc/", mode: 'copy', pattern: "*data"
  publishDir "${params.dirs.results_dir}/qc/", mode: 'copy', pattern: "*.html"

  input:
  path(aln_hc_metrics)
  path(aln_hc_dedup_metrics)

  output:
  path("*data")
  path("*.html")

  script:
  """
  multiqc --filename post_aln_filter.multiqc . 2>&1|tee > multiqc.post_aln_filter.log
  """
}


process multiqc_from_fastq_aselux {
  tag "multiqc_from_fastq_aselux"
  label "multiqc"

  publishDir "${params.dirs.results_dir}/qc/", mode: 'copy', pattern: "*data"
  publishDir "${params.dirs.results_dir}/qc/", mode: 'copy', pattern: "*.html"

  input:
  path(fastqc_raw_zip)
  path(fastqc_raw_html)
  path(trimmomatic_log)
  path(fastqc_trimmed_zip)
  path(fastqc_trimmed_html)

  output:
  path("*data")
  path("*html")

  script:
  """
  mkdir pre_trim
  mv ${fastqc_raw_zip} pre_trim
  mv ${fastqc_raw_html} pre_trim
  multiqc --filename pre_trim.multiqc ./pre_trim 2>&1|tee > multiqc.pre_trim.log
  multiqc --filename post_trim.multiqc --ignore pre_trim . 2>&1|tee > multiqc.post_trim.log
  """
}

