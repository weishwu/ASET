process exon_merge_by_gene {
    tag "exon_merge_by_gene"
    label "bedtools"

    publishDir "${params.dirs.results_dir}/ref_data/", mode: 'copy', pattern: "*bed"

    input:
    path(ref_exons_bed)
    
    output:
    path("gtf.exonMergedByGenes.bed"), emit: exons_mergedbygene_bed

    script:
    """
    awk '{OFS="\\t"; print \$1"__"\$7"__"\$6,\$2,\$3}' ${ref_exons_bed} | sortBed -i - | mergeBed -i - | sed 's/__/\\t/g' | awk '{OFS="\\t"; print \$1,\$4,\$5,\$1":"\$4+1":"\$5":"\$3":"\$2,".",\$3}' | sort -k1,1V -k2,2n >gtf.exonMergedByGenes.bed \
    2>&1|tee > exon_mergeByGene.log
    """ 
}


