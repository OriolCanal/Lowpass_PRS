process MULTIQC {
    label 'process_low'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    
    // Publish the final report to your results directory
    publishDir "${params.outdir}/00_multiqc_report", mode: 'copy'

    input:
    path 'fastp_logs/*'
    path 'picard_logs/*'

    output:
    path "*multiqc_report.html", emit: report
    path "*multiqc_data"       , emit: data

    script:
    """
    multiqc .
    """
}