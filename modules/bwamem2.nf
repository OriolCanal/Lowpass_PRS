process BWAMEM2_ALIGN {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/sangerpathogens/bwa-mem2-samtools:1.21'

    input:
    tuple val(meta), path(reads)
    path fasta
    path index_files

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
    
    script:
    def index_prefix = index_files[0].name.substring(0, index_files[0].name.lastIndexOf('.'))    
    
    """
    bwa-mem2 mem \\
        -t $task.cpus \\
        -R '@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:lib1' \\
        ${index_prefix} \\
        ${reads[0]} \\
        ${reads[1]} | \\
    samtools sort -@ $task.cpus -o ${meta.id}.sorted.bam - 

    samtools index -@ $task.cpus ${meta.id}.sorted.bam
    """
}