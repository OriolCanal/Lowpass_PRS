process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/picard:3.4.0--hdfd78af_0'

    input:  tuple val(meta), path(bam), path(bai)
    output: tuple val(meta), path("*.md.bam"), path("*.{bai,bam.bai}"), emit: bam_and_index
    tuple val(meta), path("*.metrics.txt")                    , emit: metrics // CRITICAL FOR MULTIQC

    script:
    """
    export JAVA_TOOL_OPTIONS="-Xmx${task.memory.toGiga() - 1}g"

    picard MarkDuplicates \\
        -I $bam \\
        -O ${meta.id}.md.bam \\
        -M ${meta.id}.metrics.txt \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY SILENT
    """
}