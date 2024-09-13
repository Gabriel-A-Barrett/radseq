include { BCFTOOLS_NORM            } from '../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX     } from '../../modules/nf-core/tabix/tabix/main.nf'
include { TABIX_TABIX as TABIX_FILTERS_1 } from '../../modules/nf-core/tabix/tabix/main.nf'
include { TABIX_TABIX as TABIX_FITLERS_2 } from '../../modules/nf-core/tabix/tabix/main.nf'
include { BCFTOOLS_VIEW as FILTER1 } from '../../modules/nf-core/bcftools/view/main'
include { RADSEQ_FILTERS as RADSEQ_FILTERS_1 } from '../../modules/local/bcftools/main'
include { RADSEQ_FILTERS as RADSEQ_FILTERS_2 } from '../../modules/local/bcftools/main'

workflow VCF_BCFTOOLS_RADSEQ_FILTERS {

    take:
    vcf_tbi_fasta // [[id:],vcf,tbi,fasta]

    main:
    ch_versions = Channel.empty()
    
    ch_vcf = BCFTOOLS_NORM ( vcf_tbi_fasta ).vcf
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    TABIX ( ch_vcf )

    ch_vcf_tbi = ch_vcf.join(TABIX.out.tbi)

    def fraction_missingness = params.fraction_missingness_list.toString().split(',') as List ?: [0.01,0.1]
    def minor_allele_count = params.minor_allele_count_list.toString().split(',') as List ?: [5,10,20]

    FILTER1 ( ch_vcf_tbi, [], [], [], Channel.value(fraction_missingness), Channel.value(minor_allele_count))

    RADSEQ_FILTERS_1 ( FILTER1.out.vcf , [], 'first') // site-based filtering and collect indv
    
    RADSEQ_FILTERS_2 ( FILTER1.out.vcf, RADSEQ_FILTERS_1.out.txt, 'second') // remove indv. and re-apply site-based filtering

}