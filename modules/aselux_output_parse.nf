process aselux_output_parse {
    tag "aselux_output_parse: ${sample}"
    label "pandas"

    input:
    tuple val(sample), path(ase), path(snps)
    
    output:
    path("*.aselux.hetSNP.txt"), emit: ase_count_hetsnp_txt
    path("*.aselux.homSNP.txt"), emit: ase_count_homsnp_txt
    path("*.aselux.homRef.txt"), emit: ase_count_homref_txt

    script:
    """
    aselux_output_parse.py ${sample} ${snps} ${ase} ${sample}.aselux.hetSNP.txt ${sample}.aselux.homSNP.txt ${sample}.aselux.homRef.txt \
    2>&1|tee > aselux_output_parse.${sample}.log
    """ 
}

