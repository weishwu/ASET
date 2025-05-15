include { nmask_genome } from '../modules/nmask_genome.nf'
include { star_nmaskedgenome_index } from '../modules/star_nmaskedgenome_index.nf'
include { star_nmaskedgenome_aln } from '../modules/star_nmaskedgenome_aln.nf'

workflow star_nmask_routine {
    take:
        genome
        gtf
        samples
        trimmed_reads
    
    main:
         nmask_genome(
            genome,
            samples.map { sample, read1, read2, snps, snps_with_ref -> [ sample, snps ]})

         star_nmaskedgenome_index(
            gtf,
            nmask_genome.out)

          star_nmaskedgenome_aln(
            star_nmaskedgenome_index.out.join(trimmed_reads, by: 0))

    emit:
        aln_bam = star_nmaskedgenome_aln.out.aln_bam
        aln_log = star_nmaskedgenome_aln.out.aln_log
}



