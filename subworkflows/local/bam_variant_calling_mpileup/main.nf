//
// MPILEUP variant calling: BCFTOOLS for variantcalling, SAMTools for controlfreec input
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BCFTOOLS_MPILEUP                           } from '../../../modules/nf-core/bcftools/mpileup/main'

workflow BAM_VARIANT_CALLING_MPILEUP {
    take:
    cram_crai      // channel: [mandatory] [ meta, cram, crai,]
    fasta
    intervals
    keep_bcftools_mpileup

    main:
    versions = Channel.empty()
    
    // expand channel across bed regions for variant calling multi-threading
    ch_cram_crai_bed_fasta = cram_crai
        .combine(intervals, by: [0])
        .combine(fasta, by: [0])
        .map { meta, cram, crai, bed, fasta -> 
            [[
                id:           meta.id,
                single_end:   meta.single_end,
                interval:     bed.getName().tokenize( '.' )[0],
                ref_id:       meta.ref_id
            ],
            cram, crai, file(bed), fasta]
        }

    ch_cram_crai_bed_fasta.view()

    BCFTOOLS_MPILEUP(ch_cram_crai_bed_fasta, keep_bcftools_mpileup)

    //versions = versions.mix(BCFTOOLS_MPILEUP.out.versions)

    emit:

    versions
}