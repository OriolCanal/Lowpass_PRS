process BCFTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0'

    input:  tuple val(meta), path(bam), path(bai)
            path fasta
    output: tuple val(meta), path("*.likelihoods.bcf"), emit: bcf

    script:
    """
    # Generate Genotype Likelihoods across the target regions
    bcftools mpileup \\
        -f $fasta \\
        -I -E -Q 20 -q 20 \\
        -Ob \\
        -o ${meta.id}.likelihoods.bcf \\
        --threads $task.cpus \\
        $bam

    bcftools index --threads $task.cpus ${meta.id}.likelihoods.bcf
    """
}