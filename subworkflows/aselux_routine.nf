include { aselux_index_genome } from '../modules/aselux_index_genome.nf'
include { aselux_ase } from '../modules/aselux_ase.nf'
include { aselux_output_parse } from '../modules/aselux_output_parse.nf'
include { ase_count_concat } from '../modules/ase_count_concat.nf'

workflow aselux_routine {
    take:
        genome
        gtf
        samples
    
    main:
        aselux_index_genome(genome, gtf)

        aselux_ase(
          aselux_index_genome.out.first(),
          samples.map { sample, read1, read2, snps, snps_with_ref -> [ sample, read1, read2, snps ] })
       
        aselux_output_parse(aselux_ase.out.join(samples.map { sample, read1, read2, snps, snps_with_ref -> [ sample, snps ] }, by: 0))

    emit:
        ase_hetsnp = aselux_output_parse.out.ase_count_hetsnp_txt
        ase_homsnp = aselux_output_parse.out.ase_count_homsnp_txt
        ase_homref = aselux_output_parse.out.ase_count_homref_txt
}

