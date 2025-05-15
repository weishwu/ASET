include { gmap_index_genome } from '../modules/gmap_index_genome.nf'
include { gsnap_index_gtf } from '../modules/gsnap_index_gtf.nf'
include { gsnap_index_snps } from '../modules/gsnap_index_snps.nf'
include { gsnap_aln } from '../modules/gsnap_aln.nf'
include { gsnap_sam_to_sorted_bam } from '../modules/gsnap_sam_to_sorted_bam.nf'

workflow gsnap_routine {
    take:
        genome
        gtf
        samples
        trimmed_reads
    
    main:
        gmap_index_genome(genome) 
        
        gsnap_index_gtf(gtf)

        gsnap_index_snps(
            gmap_index_genome.out.first(),
            samples.map { sample, read1, read2, snps, snps_with_ref -> [ sample, snps ] })
        
        gsnap_aln(
            gsnap_index_gtf.out.first(),
            trimmed_reads.join(gsnap_index_snps.out, by: 0))
        
        gsnap_sam_to_sorted_bam(gsnap_aln.out.aln_sam)
    
    emit:
        aln_log = gsnap_aln.out.aln_log
        aln_bam = gsnap_sam_to_sorted_bam.out
}

