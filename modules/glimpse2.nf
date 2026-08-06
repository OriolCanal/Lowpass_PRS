process GLIMPSE2_CHUNK {
    tag "$chr"
    label 'process_low'
    container 'quay.io/biocontainers/glimpse-bio:2.0.1--ha5d29c5_3'

    input:
    val chr
    path ref_sites
    path ref_sites_csi
    path genetic_map

    output:
    tuple val(chr), path("chunks.${chr}.txt"), emit: chunks

    script:
    """
    GLIMPSE2_chunk --input ${ref_sites} \\
                   --region ${chr} \\
                   --map ${genetic_map} \\
                   --sequential \\
                   --output chunks.${chr}.txt
    """
}

process GLIMPSE2_SPLIT_REFERENCE {
    tag "${chr}_${chunk_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/glimpse-bio:2.0.1--ha5d29c5_3'

    input:
    tuple val(chr), val(chunk_id), val(irg), val(org)
    path ref_panel
    path ref_panel_csi
    path genetic_map

    output:
    tuple val(chr), val(irg), val(org), path("*.bin"), emit: bin_ref

    script:
    def prefix = "split_1000GP_${chr}_chunk_${chunk_id}"
    """
    GLIMPSE2_split_reference --reference ${ref_panel} \\
                             --map ${genetic_map} \\
                             --input-region "${irg}" \\
                             --output-region "${org}" \\
                             --output ${prefix}
    """
}

process GLIMPSE2_PHASE {
    tag "${meta.id}_${chr}_${irg.replaceAll(':', '_')}"
    label 'process_high'
    container 'quay.io/biocontainers/glimpse-bio:2.0.1--ha5d29c5_3'

    input:
    tuple val(meta), path(bam), path(bai), val(chr), val(irg), val(org), path(bin_ref)

    output:
    tuple val(meta), val(chr), path("*.bcf"), path("*.bcf.csi"), emit: chunk_bcf

    script:
    def prefix = "${meta.id}_${chr}_phased_${org.replaceAll(/[:\-]/, '_')}"
    """
    GLIMPSE2_phase --bam-file ${bam} \\
                   --reference ${bin_ref} \\
                   --output ${prefix}.bcf \\
                   --threads ${task.cpus}
    """
}

process GLIMPSE2_LIGATE {
    tag "${meta.id}_${chr}"
    label 'process_medium'
    container 'quay.io/biocontainers/glimpse-bio:2.0.1--ha5d29c5_3'

    input:
    tuple val(meta), val(chr), path(chunk_bcfs), path(chunk_csis)

    output:
    tuple val(meta), path("${meta.id}_${chr}_imputed.vcf.gz"), path("${meta.id}_${chr}_imputed.vcf.gz.csi"), emit: imputed_vcf

    script:
    def list_file = "chunks_list_${meta.id}_${chr}.txt"
    def output_name = "${meta.id}_${chr}_imputed"
    """
    ls -1v ${chunk_bcfs} > ${list_file}

    GLIMPSE2_ligate --input ${list_file} \\
                    --output ${output_name}.vcf.gz \\
                    --threads ${task.cpus}
    """
}