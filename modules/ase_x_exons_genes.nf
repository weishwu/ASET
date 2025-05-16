process ase_x_exons_genes {
    tag "ase_x_exons_genes"
    label "bedtools"

    input:
    path(genes_bed)
    path(ase_hetsnp)
    path(exons_mergedbygene_bed)

    output:
    path("ASE_hetSNP_xGenes.bed"), emit: ase_x_genes_bed
    path("ASE_hetSNP_xExons.bed"), emit: ase_x_exons_bed

    script:
    """
    sed '1d' ${ase_hetsnp} | awk '{OFS="\\t";print \$1,\$2-1,\$2,".",".",\$(NF-1)}'|sed 's/plus/\\+/g' | sed 's/minus/\\-/g' | sort -k1,1V -k2,2n -k6,6V | uniq >ase_hetsnps.uniq.bed

    if [[ ${params.tool_parameters.mapper} != "ASElux" ]]; then 
       intersectBed -wa -wb -s -a ase_hetsnps.uniq.bed -b ${genes_bed} | cut -f 1,3,6,10 > ASE_hetSNP_xGenes.bed
       intersectBed -wa -wb -s -a ase_hetsnps.uniq.bed -b ${exons_mergedbygene_bed} | cut -f 1,3,6,10 > ASE_hetSNP_xExons.bed

    else 
       intersectBed -wa -wb -a ase_hetsnps.uniq.bed -b ${genes_bed} | cut -f 1,3,6,10 > ASE_hetSNP_xGenes.bed
       intersectBed -wa -wb -a ase_hetsnps.uniq.bed -b ${exons_mergedbygene_bed} | cut -f 1,3,6,10 > ASE_hetSNP_xExons.bed
    fi

    """ 
}

