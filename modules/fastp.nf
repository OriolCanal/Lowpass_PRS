process FASTP {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/fastp:1.3.3--h43da1c4_0'

    input:  tuple val(meta), path(reads)
    output: tuple val(meta), path("*.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("*.json")               , emit: json // CRITICAL FOR MULTIQC
    tuple val(meta), path("*.html")               , emit: html

    script:
    """
    fastp \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${meta.id}_R1.trimmed.fastq.gz \\
        -O ${meta.id}_R2.trimmed.fastq.gz \\
        --json ${meta.id}_fastp.json \\
        --html ${meta.id}_fastp.html \\
        --thread $task.cpus
    """
}