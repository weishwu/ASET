process trimmomatic {
    tag "trimmomatic on raw reads: $sample"
    label "trimmomatic"

    publishDir "${params.dirs.results_dir}/qc/trimmed_reads", mode: 'copy', pattern: "*.gz"
    publishDir "${params.dirs.results_dir}/qc/trimmed_reads/", mode: 'copy', pattern: "*.log"

    input: 
    tuple val(sample), path(read1), path(read2)
    path(adapters)

    output:
    tuple val(sample), path("Sample_${sample}_R1_trimmed.fastq.gz"), path("Sample_${sample}_R2_trimmed.fastq.gz"), emit: trimmomatic_reads
    path("Sample_${sample}_trimmomatic.log"), emit: trimmomatic_log

    script:
    """
    export JAVA_OPTIONS=-Xmx${task.memory.toGiga()}g
    trimmomatic PE \
        -threads ${task.cpus} \
        -${params.tool_parameters.trimmomatic.phred_version} \
        ${read1} \
        ${read2} \
        Sample_${sample}_R1_trimmed.fastq.gz \
        Sample_${sample}_R1_unpaired.fastq.gz \
        Sample_${sample}_R2_trimmed.fastq.gz \
        Sample_${sample}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:${adapters}:${params.tool_parameters.trimmomatic.ILLUMINACLIP} \
        ${params.tool_parameters.trimmomatic.more_flags} \
    >Sample_${sample}_trimmomatic.log 2>&1
    """ 
}



