include { star_genome_index } from '../modules/star_genome_index.nf'
include { star_wasp_aln } from '../modules/star_wasp_aln.nf'

workflow star_wasp_routine {
    take:
        genome
        gtf
        samples
        trimmed_reads
    
    main:
        star_genome_index(genome, gtf)

        star_wasp_aln(
          star_genome_index.out.star_genome.first(),
          trimmed_reads.join(samples.map { sample, read1, read2, snps, snps_with_ref -> [ sample, snps ] }, by: 0))
    
    emit:
        aln_bam = star_wasp_aln.out.aln_bam
        aln_log = star_wasp_aln.out.aln_log
}

