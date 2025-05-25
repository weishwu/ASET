#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

if (params.data.routine == 'from_fastq') {

 Channel
    .fromPath(params.data.sample_sheet, checkIfExists:true)
    .splitCsv(header: true)
    .map { row -> 
           def sample = row.sample
           def read1 = row.read1
           def read2 = row.read2
           def snps = row.snps
           def snps_with_ref = row.snps_with_ref ?: row.snps
           tuple(sample, file(read1), file(read2), file(snps), file(snps_with_ref))}
    .set { samples_ch }

 snps_ch = samples_ch.map { sample, read1, read2, snps, snps_with_ref -> [sample, snps_with_ref] }

} else if (params.data.routine == 'from_bam') {

 Channel
    .fromPath(params.data.sample_sheet, checkIfExists:true)
    .splitCsv(header: true)
    .map { row -> 
           def sample = row.sample
           def bam = row.bam
           def snps = row.snps
           def snps_with_ref = row.snps_with_ref ?: row.snps
           tuple(sample, file(bam), file(snps), file(snps_with_ref))}
    .set { samples_ch }
 
 snps_ch = samples_ch.map { sample, bam, snps, snps_with_ref -> [sample, snps_with_ref] }
}

genome_fa_ch = Channel.fromPath(params.support_files.genome_fa, checkIfExists:true)
gtf_ch = Channel.fromPath(params.support_files.gtf, checkIfExists:true)

if ((!params.data.routine == 'from_bam') || (params.tool_parameters.mapper != 'ASElux')) {
rRNA_tRNA_ch = Channel.fromPath(params.support_files.rRNA_tRNA_bed, checkIfExists:true)
}

include { genome_link } from './modules/genome_prep.nf'
include { genome_dict } from './modules/genome_prep.nf'
include { rRNA_tRNA_list } from './modules/rRNA_tRNA_list.nf'
include { gtf_to_df } from './modules/gtf_to_df.nf'
include { gtf_to_refflat } from './modules/gtf_to_refflat.nf'
include { extract_exons_genes } from './modules/extract_exons_genes.nf'
include { exon_merge_by_gene } from './modules/exon_merge_by_gene.nf'
include { fastqc_raw } from './modules/fastqc.nf'
include { trimmomatic } from './modules/trimmomatic.nf'
include { fastqc_trimmed } from './modules/fastqc.nf'
include { multiqc_from_fastq } from './modules/multiqc.nf'
include { multiqc_from_bam } from './modules/multiqc.nf'
include { multiqc_from_fastq_aselux } from './modules/multiqc.nf'
include { aln_filter } from './modules/aln_filter.nf'
include { aln_hc_metrics } from './modules/aln_hc_metrics.nf'
include { aln_dedup } from './modules/aln_dedup.nf'
include { aln_strand_split } from './modules/aln_strand_split.nf'
include { aln_ase_count } from './modules/aln_ase_count.nf'
include { ase_count_split } from './modules/ase_count_split.nf'
include { ase_homref_add_matGT } from './modules/ase_homref_add_matGT.nf'
include { ase_count_concat } from './modules/ase_count_concat.nf'
include { ase_x_exons_genes } from './modules/ase_x_exons_genes.nf'
include { ase_exon_label } from './modules/ase_exon_label.nf'
include { ase_hetsnp_ase } from './modules/ase_hetsnp_ase.nf'
include { ase_hetsnp_maternalcontam } from './modules/ase_hetsnp_maternalcontam.nf'
include { ase_hetsnp_add_phasing } from './modules/ase_hetsnp_add_phasing.nf'
include { ase_hetsnp_add_ase } from './modules/ase_hetsnp_add_ase.nf'
include { ase_hetsnp_add_genetype } from './modules/ase_hetsnp_add_genetype.nf'
include { qc_summary_from_fastq } from './modules/qc_summary.nf'
include { qc_summary_from_bam } from './modules/qc_summary.nf'
include { qc_summary_from_fastq_aselux } from './modules/qc_summary.nf'
include { prep_ase_data_for_plot } from './modules/prep_ase_data_for_plot.nf'
include { ase_data_rds } from './modules/ase_data_rds.nf'
include { ase_homref_x_exons_genes } from './modules/ase_homref_x_exons_genes.nf'
include { get_fastqc_version } from './modules/get_versions.nf'
include { get_ngsutils_version } from './modules/get_versions.nf'
include { get_bedtools_version } from './modules/get_versions.nf'
include { get_gatk_version } from './modules/get_versions.nf'
include { get_gsnap_version } from './modules/get_versions.nf'
include { get_multiqc_version } from './modules/get_versions.nf'
include { get_star_version } from './modules/get_versions.nf'
include { get_trimmomatic_version } from './modules/get_versions.nf'
include { get_aselux_version } from './modules/get_versions.nf'
include { get_pandas_version } from './modules/get_versions.nf'
include { get_r_version } from './modules/get_versions.nf'
include { version_report } from './modules/version_report.nf'

include { gsnap_routine } from './subworkflows/gsnap_routine.nf'
include { star_wasp_routine } from './subworkflows/star_wasp_routine.nf'
include { star_nmask_routine } from './subworkflows/star_nmask_routine.nf'
include { aselux_routine } from './subworkflows/aselux_routine.nf'


workflow {

    genome_fa_ch | genome_link | genome_dict
    gtf_ch | gtf_to_df | extract_exons_genes  
    gtf_ch | gtf_to_refflat
    extract_exons_genes.out.exons_bed | exon_merge_by_gene

 if (params.data.routine == 'from_fastq') {

    fastqc_raw(samples_ch.map { sample, read1, read2, snps, snps_with_ref -> [ sample, read1, read2 ] })    
    
    adapters = channel.fromPath(params.support_files.adapters, checkIfExists:true)
    trimmomatic(
       samples_ch.map { sample, read1, read2, snps, snps_with_ref -> [ sample, read1, read2 ] },
       adapters.first())

    fastqc_trimmed(trimmomatic.out.trimmomatic_reads)
 }

 if ((params.tool_parameters.mapper == "STAR_WASP") && (params.data.routine == "from_fastq")) {

    star_wasp_routine(
       genome_link.out,
       gtf_ch,
       samples_ch,
       trimmomatic.out.trimmomatic_reads)

    raw_aln = star_wasp_routine.out.aln_bam
    aln_log = star_wasp_routine.out.aln_log
  }
  
 if ((params.tool_parameters.mapper == "GSNAP") && (params.data.routine == "from_fastq")) {

    gsnap_routine(
       genome_link.out,
       gtf_ch,
       samples_ch,
       trimmomatic.out.trimmomatic_reads)

    raw_aln = gsnap_routine.out.aln_bam
    aln_log = gsnap_routine.out.aln_log
  }
  
 if ((params.tool_parameters.mapper == "STAR_NMASK") && (params.data.routine == "from_fastq")) {
  
    star_nmask_routine(
       genome_link.out.first(),
       gtf_ch.first(),
       samples_ch,
       trimmomatic.out.trimmomatic_reads)

    raw_aln = star_nmask_routine.out.aln_bam
    aln_log = star_nmask_routine.out.aln_log 
  }
 
 if (params.data.routine == 'from_bam') {
    raw_aln = samples_ch.map { sample, bam, snps, snps_with_ref -> [sample, bam] }
 }

 if (params.tool_parameters.mapper != "ASElux") {
    aln_filter(raw_aln)

    rRNA_tRNA_list(rRNA_tRNA_ch, genome_dict.out)

    aln_hc_metrics(
       aln_filter.out.aln_hc_bam,
       gtf_to_refflat.out.first(),
       rRNA_tRNA_list.out.first())

    aln_dedup(aln_filter.out.aln_hc_bam)
 
    aln_strand_split(aln_dedup.out.aln_hc_dedup_bam)

    aln_ase_count(
       genome_link.out.first(), 
       genome_dict.out.first(), 
       aln_strand_split.out.aln_strand_split_bam.join(snps_ch, by: 0))
    
       ase_count_split(
       aln_ase_count.out.join(snps_ch, by: 0))

    ase_split_hetsnp = ase_count_split.out.ase_count_hetsnp_txt.collect()
    ase_split_homsnp = ase_count_split.out.ase_count_homsnp_txt.collect()
    ase_split_homref = ase_count_split.out.ase_count_homref_txt.collect()

 } else {
       aselux_routine(
       genome_link.out.first(),
       gtf_ch.first(),
       samples_ch,
       trimmomatic.out.trimmomatic_reads)

    ase_split_hetsnp = aselux_routine.out.ase_hetsnp.collect()
    ase_split_homsnp = aselux_routine.out.ase_homsnp.collect()
    ase_split_homref = aselux_routine.out.ase_homref.collect()
 }

    ase_count_concat(
       ase_split_hetsnp, 
       ase_split_homsnp, 
       ase_split_homref)

    ase_hetsnp_ase(
       ase_count_concat.out.ase_count_hetsnp_txt,
       ase_count_concat.out.ase_count_homsnp_txt,
       ase_count_concat.out.ase_count_homref_txt)

    ase_x_exons_genes(
       extract_exons_genes.out.genes_bed,
       ase_count_concat.out.ase_count_hetsnp_txt,
       exon_merge_by_gene.out.exons_mergedbygene_bed)

    ase_exon_label(
       ase_hetsnp_ase.out,
       ase_x_exons_genes.out.ase_x_exons_bed,
       ase_x_exons_genes.out.ase_x_genes_bed)

    ase_hetsnp_add_genetype(
       ase_exon_label.out,
       extract_exons_genes.out.genes_bed)

    ase_hetsnp = ase_hetsnp_add_genetype.out

// if there is phasing data, use it to assign parental allelic counts
  
  if (params.data.phase_info) {

    Channel
      .fromPath(params.data.phase_info, checkIfExists:true)
      .splitCsv(header: true)
      .filter { row -> row.phase_info != '.' }
      .map { row ->
          def sample = row.sample
          def phase_info = row.phase_info
          tuple(sample, phase_info)}
      .set { phase_info }

    phase_info_files = phase_info.map { sample, phase_info -> [phase_info]}
    ase_hetsnp_add_phasing(
       ase_hetsnp,
       Channel.fromPath(params.data.phase_info, checkIfExists:true),
       phase_info_files.collect().flatten().unique())

    ase_hetsnp = ase_hetsnp_add_phasing.out
  }

// if there is maternal_vcf, add gene-level maternal contamination data (for placenta samples)    
  if (params.data.maternal_vcf) {

    Channel
       .fromPath(params.data.maternal_vcf, checkIfExists:true)
       .splitCsv(header: true)
       .filter { row -> row.mat_snps != '.' }
       .map { row ->
              def sample = row.sample
              def mat_snps = row.mat_snps
              tuple(sample, mat_snps)}
       .set { mat_vcf }

    mat_vcf_files = mat_vcf.map { sample, mat_snps -> [mat_snps]}
   
    ase_homref_add_matGT(
       ase_count_concat.out.ase_count_homref_txt,
       Channel.fromPath(params.data.maternal_vcf, checkIfExists:true),
       mat_vcf_files.collect().flatten().unique())

    ase_homref_x_exons_genes(
       ase_homref_add_matGT.out,
       exon_merge_by_gene.out.exons_mergedbygene_bed)

    ase_hetsnp_maternalcontam(
       ase_hetsnp,
       ase_homref_x_exons_genes.out)

    ase_hetsnp = ase_hetsnp_maternalcontam.out
  }

    ase_data_rds(
       ase_hetsnp, 
       exon_merge_by_gene.out.exons_mergedbygene_bed)

 if ((params.data.routine == 'from_fastq') && (params.tool_parameters.mapper != "ASElux")) {

    multiqc_from_fastq(
           fastqc_raw.out.fastqc_zip.collect(),
           fastqc_raw.out.fastqc_html.collect(),
           trimmomatic.out.trimmomatic_log.collect(),
           fastqc_trimmed.out.fastqc_zip.collect(),
           fastqc_trimmed.out.fastqc_html.collect(),
           aln_log.collect(),
           aln_hc_metrics.out.collect(),
           aln_dedup.out.aln_hc_dedup_metrics.collect())

    qc_summary_from_fastq(
           ase_hetsnp,
           Channel.fromPath(params.data.sample_sheet, checkIfExists:true),
           trimmomatic.out.trimmomatic_log.collect(),
           aln_log.collect(),
           aln_filter.out.aln_hc_flagstat.collect(),
           aln_strand_split.out.aln_hc_dedup_flagstat.collect(),
           aln_hc_metrics.out.collect(),
           samples_ch.map {sample, read1, read2, snps, snps_with_ref -> [snps]}.collect())
 }


 if ((params.data.routine == 'from_fastq') && (params.tool_parameters.mapper == "ASElux")) {

    multiqc_from_fastq_aselux(
           fastqc_raw.out.fastqc_zip.collect(),
           fastqc_raw.out.fastqc_html.collect(),
           trimmomatic.out.trimmomatic_log.collect(),
           fastqc_trimmed.out.fastqc_zip.collect(),
           fastqc_trimmed.out.fastqc_html.collect())

    qc_summary_from_fastq_aselux(
           ase_hetsnp,
           Channel.fromPath(params.data.sample_sheet, checkIfExists:true),
           trimmomatic.out.trimmomatic_log.collect(),
           samples_ch.map {sample, read1, read2, snps, snps_with_ref -> [snps]}.collect())
 }

 if (params.data.routine == 'from_bam') {
 
    multiqc_from_bam(
           aln_hc_metrics.out.collect(),
           aln_dedup.out.aln_hc_dedup_metrics.collect())

    qc_summary_from_bam(
           ase_hetsnp,
           Channel.fromPath(params.data.sample_sheet, checkIfExists:true),
           aln_filter.out.aln_hc_flagstat.collect(),
           aln_strand_split.out.aln_hc_dedup_flagstat.collect(),
           aln_hc_metrics.out.collect(),
           samples_ch.map {sample, bam, snps, snps_with_ref -> [snps]}.collect())
 }

    get_aselux_version()
    get_bedtools_version()
    get_fastqc_version()
    get_gatk_version()
    get_gsnap_version()
    get_multiqc_version()
    get_ngsutils_version()
    get_pandas_version()
    get_r_version()
    get_star_version()
    get_trimmomatic_version()
    
 if ((params.data.routine == 'from_fastq') && ((params.tool_parameters.mapper == "STAR_WASP") || (params.tool_parameters.mapper == "STAR_NMASK"))) {
    version_report(
        get_star_version.out.concat(
        get_gatk_version.out).concat(
        get_fastqc_version.out).concat(
        get_trimmomatic_version.out).concat(
        get_bedtools_version.out).concat(
        get_multiqc_version.out).concat(
        get_pandas_version.out).concat(
        get_r_version.out).collect())
 }

 if ((params.data.routine == 'from_fastq') && (params.tool_parameters.mapper == "GSNAP")) {
    version_report(
        get_gsnap_version.out.concat(
        get_gatk_version.out).concat(
        get_fastqc_version.out).concat(
        get_trimmomatic_version.out).concat(
        get_bedtools_version.out).concat(
        get_multiqc_version.out).concat(
        get_pandas_version.out).concat(
        get_r_version.out).collect())
 }

 if ((params.data.routine == 'from_fastq') && (params.tool_parameters.mapper == "ASElux")) {
    version_report(
        get_fastqc_version.out.concat(
        get_trimmomatic_version.out).concat(
        get_bedtools_version.out).concat(
        get_multiqc_version.out).concat(
        get_pandas_version.out).concat(
        get_r_version.out).collect())
 }

 if (params.data.routine == 'from_bam') {
    multiqc_version_report(
        get_bedtools_version.out.concat(
        get_multiqc_version.out).concat(
        get_pandas_version.out).concat(
        get_r_version.out).collect())
 }

}

