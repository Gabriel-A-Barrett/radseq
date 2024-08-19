//
// MPILEUP variant calling: BCFTOOLS for variantcalling
//

include { BCFTOOLS_MPILEUP                           } from '../../../modules/nf-core/bcftools/mpileup/main'
include { TABIX_TABIX } from '../../../modules/nf-core/tabix/tabix/main'
include { VCF_GATHER_BCFTOOLS } from '../../nf-core/vcf_gather_bcftools/main'

workflow BAM_VARIANT_CALLING_MPILEUP {
    take:
    cram_crai      // channel: [mandatory] [ meta, cram, crai,]
    fasta          // channel: [mandatory] [ meta, fasta ]
    intervals      // channel: [mandatory] [ meta, intervals]
    keep_bcftools_mpileup // boolean: save intermediate .mpileup files

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
            cram, crai, file(bed), fasta ]
        }

    BCFTOOLS_MPILEUP(ch_cram_crai_bed_fasta, keep_bcftools_mpileup)
    versions = versions.mix(BCFTOOLS_MPILEUP.out.versions)

    //TABIX_TABIX ( BCFTOOLS_MPILEUP.out.bcf )

    // channel
    ch_bcf_tbi = BCFTOOLS_MPILEUP.out.bcf.combine( BCFTOOLS_MPILEUP.out.csi, by: [0] )

    ch_bcf = BCFTOOLS_MPILEUP.out.bcf
    ch_csi = BCFTOOLS_MPILEUP.out.csi

    // is there more than one output channel
    if (ch_bcf_tbi.count().map { it > 1 } ) { 
        
        ch_mpile_bcfs = ch_bcf_tbi
            .map{ meta, bcf, tbi -> 
            [[
                id: meta.id,
                ref_id: meta.ref_id
            ], 
            bcf, tbi ] 
            }.groupTuple()
        
        ch_mpile_bcfs.view()
        //
        // SUBWORKFLOW
        //
        VCF_GATHER_BCFTOOLS ( ch_mpile_bcfs, [], true)
        
        // Overwrite channels
        ch_bcf = VCF_GATHER_BCFTOOLS.out.vcf
        ch_csi = VCF_GATHER_BCFTOOLS.out.csi
    }

    emit:
    bcf = ch_bcf
    tbi = ch_csi

    versions
}