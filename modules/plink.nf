process PLINK2_PRS {
    tag "$meta.id"
    label 'process_low'
    container 'biocontainers/plink2:v2.00a5_cv1'
    publishDir "${params.outdir}/prs_scores", mode: 'copy'

    input:
    tuple val(meta), path(vcf)
    path prs_weights

    output:
    path "*.profile", emit: scores

    script:
    """
    # 1. Filter out low-confidence imputed variants (INFO/DR2 score < 0.8)
    # 2. Calculate the calculated PRS score matrix
    plink2 --vcf $vcf dosage=DS \\
           --extract-if-info "INFO >= 0.8" \\
           --score $prs_weights 1 2 3 header \\
           --out ${meta.id}_prs
    """
}