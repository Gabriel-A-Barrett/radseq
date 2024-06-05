//
// MPILEUP variant calling: BCFTOOLS for variantcalling, SAMTools for controlfreec input
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BCFTOOLS_MPILEUP                           } from '../../../modules/nf-core/bcftools/mpileup/main'

workflow BAM_VARIANT_CALLING_MPILEUP {
    take:
    cram_crai_fasta      // channel: [mandatory] [ meta, cram, crai,]
    keep_bcftools_mpileup

    main:
    versions = Channel.empty()

    
    
    BCFTOOLS_MPILEUP(cram_crai_fasta, keep_bcftools_mpileup)

    versions = versions.mix(BCFTOOLS_MPILEUP.out.versions)

    emit:

    versions
}