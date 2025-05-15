process gsnap_index_gtf {
    tag "gsnap_index_gtf"
    label "GSNAP"

    input:
    path(gtf)

    output:
    path("trx_splicesites.iit"), emit: "gsnps_trx_index"

    script:
    """
    cat ${gtf} | gtf_splicesites > trx.splicesites
    cat trx.splicesites | iit_store -o trx_splicesites
    """
}

