process BWAMEM2_MEM {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    input:
    tuple val(meta), path(reads), path(index), path(fasta)
    val   sequence_type
    val   lengths

    output:
    tuple val(meta), path("*.bam")  , emit: bam,    optional: true
    tuple val(meta), path("*.cram") , emit: cram,   optional: true
    tuple val(meta), path("*.csi")  , emit: csi,    optional: true
    tuple val(meta), path("*.crai") , emit: crai,   optional: true
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args3.contains("--output-fmt sam")   ? "sam" :
                    args3.contains("--output-fmt cram")  ? "cram":
                    args3.contains("-O cram")            ? "cram":
                    args3.contains("-C")                 ? "cram":
                    "bam"
    def reference = fasta && extension=="cram"  ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"
    if (sequence_type == 'PE' && params.method == 'denovo') {  
        """
        INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
        
        echo "${lengths.join("\n")}" > lengths.txt
        MLEN=\$(awk '{ print length() | "sort -rn" }' lengths.txt | head -1)
        INSERT=\$(( \$MLEN * 2 ))
        INSERTH=\$(( \$INSERT + 100 ))
        INSERTL=\$(( \$INSERT - 100 ))
        SD=\$(( \$INSERT / 5 ))

        bwa-mem2 \\
            mem \\
            $args \\
            -I \$INSERT,\$SD,\$INSERTH,\$INSERTL \\
            -R "@RG\\tID:${prefix}\\tSM:${prefix}\\tPL:Illumina" \\
            -t $task.cpus \\
            \$INDEX \\
            $reads \\
            | samtools view $args2 \\
            | samtools sort $args3 $reference --threads $task.cpus -o ${prefix}.${extension} -

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS
        """
    } else {
        """
        INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

        bwa-mem2 \\
            mem \\
            $args \\
            -R "@RG\\tID:${prefix}\\tSM:${prefix}\\tPL:Illumina" \\
            -t $task.cpus \\
            \$INDEX \\
            $reads \\
            | samtools view $args2 \\
            | samtools sort $args3 $reference --threads $task.cpus -o ${prefix}.${extension} -

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS
        """        
    }
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}