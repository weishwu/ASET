process version_report {
    tag "version_report"

    publishDir "${params.dirs.results_dir}/", mode: 'copy', pattern: "tool_versions.txt"

    input:
    path(version_files)

    output:
    path("tool_versions.txt")

    script:
    """
    cat *.txt > tool_versions.txt
    """
}

