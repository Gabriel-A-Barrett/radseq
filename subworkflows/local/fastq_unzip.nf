include { GUNZIP as GUNZIP_FASTA            } from '../../modules/nf-core/gunzip/main'

workflow FASTQ_UNZIP {
    take:
    fasta

    main:
    
    ch_versions = Channel.empty()

    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta    = Channel.fromPath(fasta)
    }

    emit:
    fasta = ch_fasta

}