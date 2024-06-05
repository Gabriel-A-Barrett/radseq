//
// Alignment with BWA
//

include { BWA_MEM                 } from '../../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX          } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS      } from '../bam_stats_samtools/main'

workflow FASTQ_ALIGN_BWA {
    take:
    ch_reads_index_fasta        // channel (mandatory): [ val(meta), [ path(reads) ] ]
    fasta        // channel (optional) : [ val(meta3), path(fasta) ]
    sequence_type
    read_lengths

    main:
    ch_versions = Channel.empty()

    //
    // Map reads with BWA, filter, sort and write to cram
    //

    BWA_MEM ( ch_reads_index_fasta, sequence_type, read_lengths )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    SAMTOOLS_INDEX ( BWA_MEM.out.cram )

    //
    // run samtools stats, flagstat and idxstats
    //

    ch_cram_crai = BWA_MEM.out.cram
        .join( SAMTOOLS_INDEX.out.crai, by: [0], remainder: true )

    
    ch_cram_crai_fasta = ch_cram_crai
        .map { meta, cram, crai ->
            [[id:meta.id.split(/[^\p{L}]/)[0], single_end:meta.single_end, ref_id:meta.ref_id], cram, crai] 
        }
        .combine( fasta, by: [0] )
        .map { meta, cram, crai, fasta ->
            [[id:cram.getSimpleName(), single_end:meta.single_end, ref_id:meta.ref_id], cram, crai, fasta]
            }

    BAM_STATS_SAMTOOLS ( ch_cram_crai_fasta )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    cram      = BWA_MEM.out.cram
    cram_crai = ch_cram_crai
    cram_crai_fasta = ch_cram_crai_fasta
    flagstat  = BAM_STATS_SAMTOOLS.out.flagstat
    //bam      = BAM_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), path(bam) ]
    //bai      = BAM_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), path(bai) ]
    //csi      = BAM_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), path(csi) ]
    stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]
    versions = ch_versions                          // channel: [ path(versions.yml) ]
}
