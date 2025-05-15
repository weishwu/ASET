process ase_homref_x_exons_genes {
    tag "ase_homref_x_exons_genes"
    label "bedtools"

    input:
    path(ase_homref)
    path(exons_mergedbygene_bed)
    
    output:
    path("ASE_homRef_MatGT_exonic.gz"), emit: ase_homref_x_exons

    script:
    """
    echo -e "gene\\tnonRefFreq\\tcontig\\tposition\\tRNAid\\tMatGT" >ASE_homRef_MatGT_exonic
    ase_homref_lc=`gunzip -c ${ase_homref} | sed '1d' | wc -l`
    if [[ \${ase_homref_lc} -gt 0 ]]; then
    gunzip -c ${ase_homref} | sed '1d' | awk -v t=${params.data.ase_depth_min} '{OFS="\\t";if ((\$8+\$12) >= t) {print \$1,\$2-1,\$2,".",(1-\$6/(\$8+\$12)),\$14,\$15,\$NF}}' | sed 's/plus/\\+/g' | sed 's/minus/\\-/g' | intersectBed -wa -wb -s -a - -b ${exons_mergedbygene_bed} | awk '{OFS="\\t";print \$12,\$5,\$1,\$3,\$7,\$8}' | sed 's/:/\\t/g' | cut -f5,7- | sed 's/ //g' | sed 's:|:\\/:g' >>ASE_homRef_MatGT_exonic
    fi
    gzip ASE_homRef_MatGT_exonic
    """ 
}

