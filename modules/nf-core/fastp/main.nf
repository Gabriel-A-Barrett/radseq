process FASTP {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::fastp=0.24.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.24.0--heae3180_1' :
        'quay.io/biocontainers/fastp:0.24.0--heae3180_1' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(uniq_full_fasta)
    val   save_trimmed_fail
    val   save_merged
    val   denovo_construction

    output:
    tuple val(meta), path('*.fastp.fastq.gz') , optional:true, emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.log')            , emit: log
    path "versions.yml"                       , emit: versions
    tuple val(meta), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail
    tuple val(meta), path('*.merged.fastq.gz'), optional:true, emit: reads_merged
    tuple val(meta), path('*.uniq.fasta')     , optional:true, emit: fasta
    tuple val(meta), path('*.totaluniqseq')   , optional:true, emit: totaluniqseq
    //stdout emit: tosystem

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fail_fastq = save_trimmed_fail && meta.single_end ? "--failed_out ${prefix}.fail.fastq.gz" : save_trimmed_fail && !meta.single_end ? "--unpaired1 ${prefix}_1.fail.fastq.gz --unpaired2 ${prefix}_2.fail.fastq.gz" : ''
    // Added soft-links to original fastqs for consistent naming in MultiQC
    // Use single ended for interleaved. Add --interleaved_in in config.
    if ( task.ext.args?.contains('--interleaved_in') ) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz
        fastp \\
            --stdout \\
            --in1 ${prefix}.fastq.gz \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $fail_fastq \\
            $args \\
            2> ${prefix}.fastp.log \\
        | gzip -c > ${prefix}.fastp.fastq.gz
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else if (denovo_construction && meta.single_end) {
        """
        MaxLen="\$(awk '!/>/' ${uniq_full_fasta}  | \\
            awk '(NR==1||length<shortest){shortest=length} END {print shortest}')"
        
        fastp \\
            --in1 ${reads} \\
            --out1 ${prefix}.fastp.fastq.gz \\
            --thread ${task.cpus} \\
            ${args} \\
            --disable_quality_filtering \\
            --length_required \$MaxLen \\
            &> fastp.log
        
        # Fastq back to Fasta
        gunzip ${prefix}.fastp.fastq.gz
        awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' ${prefix}.fastp.fastq | \\
            paste - - | \\
            sort -k1,1 -V | \\
            tr "\\t" "\\n" > ${prefix}.uniq.fasta
        
        awk '!/>/' ${prefix}.uniq.fasta > ${prefix}.totaluniqseq
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
            BusyBox: \$(busybox | sed -n -E 's/.*v([[:digit:].]+)\\s\\(.*/\\1/p')
        END_VERSIONS
        """
    } else if (denovo_construction && !meta.single_end) {
        """
        MaxLen="\$(awk '!/>/' ${uniq_full_fasta}  | \\
            awk '(NR==1||length<shortest){shortest=length} END {print shortest}')"
        
        fastp \\
            --in1 ${reads} \\
            --out1 ${prefix}.fastp.fastq.gz \\
            --thread ${task.cpus} \\
            ${args} \\
            --disable_quality_filtering \\
            --length_required \$MaxLen \\
            &> fastp.log
        
        # Fastq back to Fasta
        gunzip ${prefix}.fastp.fastq.gz
        awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' ${prefix}.fastp.fastq | \\
            paste - - | \\
            sort -k1,1 -V | \\
            tr "\\t" "\\n" > ${prefix}.uniq.fasta
        
        awk '!/>/' ${prefix}.uniq.fasta > ${prefix}.totaluniqseq
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
            BusyBox: \$(busybox | sed -n -E 's/.*v([[:digit:].]+)\\s\\(.*/\\1/p')
        END_VERSIONS
        """    

    } else if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf ${reads} ${prefix}.fastq.gz
        fastp \\
            --stdout \\
            --in1 ${prefix}.fastq.gz \\
            --out1  ${prefix}.fastp.fastq.gz \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $fail_fastq \\
            $args \\
            2> ${prefix}.fastp.log
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else if (!meta.single_end && meta.umi_barcodes) {
        def merge_fastq = save_merged ? "-m --merged_out ${prefix}.merged.fastq.gz" : ''
        def umi_args = task.ext.umi_args ?: ''
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
        fastp \\
            --in1 ${prefix}_1.fastq.gz \\
            --in2 ${prefix}_2.fastq.gz \\
            --out1 ${prefix}_1.fastp.fastq.gz \\
            --out2 ${prefix}_2.fastp.fastq.gz \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $fail_fastq \\
            $merge_fastq \\
            --thread $task.cpus \\
            $args \\
            $umi_args \\
            2> ${prefix}.fastp.log
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else if (meta.single_end && meta.umi_barcodes) {
        def umi_args = task.ext.umi_args ?: '' 
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz
        fastp \\
            --stdout \\
            --in1 ${prefix}.fastq.gz \\
            --out1  ${prefix}.fastp.fastq.gz \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $fail_fastq \\
            $umi_args \\
            $args \\
            2> ${prefix}.fastp.log
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """  
    } else {
        def merge_fastq = save_merged ? "-m --merged_out ${prefix}.merged.fastq.gz" : ''
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
        fastp \\
            --in1 ${prefix}_1.fastq.gz \\
            --in2 ${prefix}_2.fastq.gz \\
            --out1 ${prefix}_1.fastp.fastq.gz \\
            --out2 ${prefix}_2.fastp.fastq.gz \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $fail_fastq \\
            $merge_fastq \\
            --thread $task.cpus \\
            $args \\
            2> ${prefix}.fastp.log
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    }
}

