process ase_count_concat {
    tag "ase_count_concat"

    publishDir "${params.dirs.results_dir}/ase_data/", mode: 'copy', pattern: "*txt"

    input:
    path(hetsnp_ase)
    path(homsnp_ase)
    path(homref_ase)    

    output:
    path("ASE_hetSNP.txt"), emit: ase_count_hetsnp_txt
    path("ASE_homSNP.txt"), emit: ase_count_homsnp_txt
    path("ASE_homRef.txt"), emit: ase_count_homref_txt

    script:
    """
    fn=`ls *.txt | head -1`
    head -1 \${fn} >ASE_hetSNP.txt
    head -1 \${fn} >ASE_homSNP.txt
    head -1 \${fn} >ASE_homRef.txt

    hetSNP_count=`cat ${hetsnp_ase} | grep -v '^contig' | wc -l`
    [ \${hetSNP_count} -eq 0 ] || cat ${hetsnp_ase} | grep -v '^contig' >>ASE_hetSNP.txt

    homSNP_count=`cat ${homsnp_ase} | grep -v '^contig' | wc -l`
    [ \${homSNP_count} -eq 0 ] || cat ${homsnp_ase} | grep -v '^contig' >>ASE_homSNP.txt

    homRef_count=`cat ${homref_ase} | grep -v '^contig'| wc -l`
    [ \${homRef_count} -eq 0 ] || cat ${homref_ase} | grep -v '^contig' >>ASE_homRef.txt
    """ 
}

