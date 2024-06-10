include { BCFTOOLS_CONCAT    } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT      } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX        } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_GATHER_BCFTOOLS {

    take:
    ch_vcfs_tbis        // channel: [ meta, vcf, tbi ]
    ch_scatter_output   // channel: [ meta, bed, gather_count ] => output from the scatter subworkflow, if you didn't use this subworkflow you can just use `[]` as bed since it isn't used
    val_sort            // boolean: Whether or not the output file should be sorted !! Add the config when using sort !!

    main:

    ch_versions = Channel.empty()

    BCFTOOLS_CONCAT ( ch_vcfs_tbis )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    if (val_sort) {
        BCFTOOLS_SORT(BCFTOOLS_CONCAT.out.vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

        ch_tabix_input = BCFTOOLS_SORT.out.vcf

    } else {
        ch_tabix_input = BCFTOOLS_CONCAT.out.vcf
    }

    TABIX_TABIX ( ch_tabix_input )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
    vcf      = ch_tabix_input        // channel: [ val(meta), [ vcf ] ]
    tbi      = TABIX_TABIX.out.tbi   // channel: [ val(meta), [ tbi ] ]
    csi      = TABIX_TABIX.out.csi

    versions = ch_versions           // channel: [ versions.yml ]
}

