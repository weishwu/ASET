process gtf_to_refflat {
  tag "gtf_to_refflat"

  input:
  path(gtf)

  output:
  path("gtf.refFlat")

  script:
  """
  gtfToGenePred -genePredExt -geneNameAsName2 ${gtf} gtf.genePred
  awk '{OFS="\\t";print \$12,\$0}' gtf.genePred | cut -f1-11 > gtf.refFlat
  rm gtf.genePred
  """
}

