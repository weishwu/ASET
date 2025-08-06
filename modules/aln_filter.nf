process aln_filter {
   tag "aln_filter: ${sample}"
   label "ngsutils"

   publishDir "${params.dirs.results_dir}/aln/02.aln_hc/", mode: 'copy', pattern: "*ba?"
   publishDir "${params.dirs.results_dir}/aln/02.aln_hc/", mode: 'copy', pattern: "*flagstat"

   input:
   tuple val(sample), path(aln)
    
   output:
   tuple val(sample), path("${sample}.coord_sorted.hc.bam"), emit: aln_hc_bam
   tuple val(sample), path("${sample}.coord_sorted.hc.bam.bai"), emit: aln_hc_bai
   path("${sample}.coord_sorted.hc.bam.flagstat"), emit: aln_hc_flagstat
   path("aln.hc_filt.${sample}.log")

   script:
   """
   half_cpus=`echo ${task.cpus} | awk '{print \$1/2}'`

   if [[ ${params.tool_parameters.mapper} == "STAR_WASP" ]]; then 
      samtools view -@${task.cpus} -bh ${params.tool_parameters.aln_filter_flags} ${aln} > ${sample}.coord_sorted.filt1.bam 2>aln.hc_filt.${sample}.log

      samtools view -H ${sample}.coord_sorted.filt1.bam >${sample}.coord_sorted.filt2.sam 2>>aln.hc_filt.${sample}.log
      samtools view -@${task.cpus} ${sample}.coord_sorted.filt1.bam | grep -v "vW:i:" >>${sample}.coord_sorted.filt2.sam 2>>aln.hc_filt.${sample}.log

      bamutils filter ${sample}.coord_sorted.filt1.bam ${sample}.coord_sorted.filt3.bam -eq vW 1 2>&1|tee >>aln.hc_filt.${sample}.log
      samtools view -@${task.cpus} ${sample}.coord_sorted.filt3.bam >>${sample}.coord_sorted.filt2.sam 2>>aln.hc_filt.${sample}.log

      samtools view -@\${half_cpus} -b ${sample}.coord_sorted.filt2.sam | samtools sort -@\${half_cpus} -O BAM -o ${sample}.coord_sorted.hc.bam 2>&1|tee >>aln.hc_filt.${sample}.log

   else
      samtools view -@${task.cpus} -bh ${params.tool_parameters.aln_filter_flags} ${aln} > ${sample}.coord_sorted.hc.bam 2>&1|tee >>aln.hc_filt.${sample}.log
   fi

   samtools index -@${task.cpus} ${sample}.coord_sorted.hc.bam 2>>aln.hc_filt.${sample}.log
   samtools flagstat -@${task.cpus} ${sample}.coord_sorted.hc.bam > ${sample}.coord_sorted.hc.bam.flagstat 2>>aln.hc_filt.${sample}.log

   rm -rf ${sample}.coord_sorted.filt1.bam  ${sample}.coord_sorted.filt2.sam ${sample}.coord_sorted.filt3.bam
   """ 
}

