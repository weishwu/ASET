process get_aselux_version {
    tag "ASElux version"
    label "ASElux"

    output:
    path("aselux.version.txt")

    script:
    """
    echo -e "ASElux\\t1.0.2" > aselux.version.txt
    """
}

process get_bedtools_version {
    tag "bedtools version"
    label "bedtools"

    output:
    path("bedtools.version.txt")

    script:
    """
    bedtools --version | head -1 | awk '{OFS="\\t"; print \$1,\$2}' > bedtools.version.txt
    """
}

process get_fastqc_version {
    tag "fastqc version"
    label "fastqc"

    output:
    path("fastqc.version.txt")

    script:
    """
    fastqc -v | head -1 | awk '{OFS="\\t"; print \$1,\$2}' > fastqc.version.txt
    """
}

process get_gatk_version {
    tag "gatk version"
    label "gatk"

    output:
    path("gatk.version.txt")

    script:
    """
    gatk --version | grep GATK | head -1 | awk '{OFS="\\t"; print "GATK",\$NF}' > gatk.version.txt
    """
}

process get_gsnap_version {
    tag "GSNAP version"
    label "GSNAP"

    output:
    path("gsnap.version.txt")

    script:
    """
    gsnap --version | grep "version" | head -1 | awk '{OFS="\\t"; print "GSNAP",\$NF}' > gsnap.version.txt
    """
}

process get_multiqc_version {
    tag "multiqc version"
    label "multiqc"

    output:
    path("multiqc.version.txt")

    script:
    """
    multiqc --version | head -1 | awk '{OFS="\\t"; print "MultiQC",\$NF}' > multiqc.version.txt
    """
}

process get_ngsutils_version {
    tag "ngsutils version"
    label "ngsutils"

    output:
    path("ngsutils.version.txt")

    script:
    """
    samtools --version | head -1 | awk '{OFS="\\t"; print \$1,\$2}' > ngsutils.version.txt 
    conda list | grep ngsutils | head -1 | awk '{OFS="\\t"; print \$1,\$2}' >> ngsutils.version.txt
    """
}

process get_pandas_version {
    tag "pandas version"
    label "pandas"

    output:
    path("pandas.version.txt")

    script:
    """
    python -V | head -1 | awk '{OFS="\\t"; print \$1,\$2}' > pandas.version.txt
    python -c "import pandas as pd; print(pd.__version__)" | head -1 | awk '{OFS="\\t"; print "pandas",\$NF}' >> pandas.version.txt
    """
}

process get_r_version {
    tag "R version"
    label "rplot"

    output:
    path("r.version.txt")

    script:
    """
    R --version | head -1 | awk '{OFS="\\t"; print \$1,\$3}' > r.version.txt
    """
}

process get_star_version {
    tag "STAR version"
    label "STAR"

    output:
    path("star.version.txt")

    script:
    """
    STAR --version | head -1 | awk '{OFS="\\t"; print "STAR",\$NF}' > star.version.txt
    """
}

process get_trimmomatic_version {
    tag "trimmomatic version"
    label "trimmomatic"

    output:
    path("trimmomatic.version.txt")

    script:
    """
    trimmomatic -version | head -1 | awk '{OFS="\\t"; print "Trimmomatic",\$NF}' > trimmomatic.version.txt
    """
}

