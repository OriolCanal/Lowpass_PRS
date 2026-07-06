#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Include modules
include { FASTP                  } from './modules/fastp'
include { MULTIQC                } from './modules/multiqc'
include { BWAMEM2_ALIGN          } from './modules/bwamem2'
include { PICARD_MARKDUPLICATES  } from './modules/picard'
include { BCFTOOLS_MPILEUP       } from './modules/bcftools'
include { GLIMPSE2_CHUNK; 
          GLIMPSE2_SPLIT_REFERENCE; 
          GLIMPSE2_PHASE; 
          GLIMPSE2_LIGATE        } from './modules/glimpse2'// include { PLINK2_PRS             } from './modules/plink'

workflow {
    // 1. Create input channel from your folder structure
    // Automatically extracts sample ID from filename or path
    ch_input = Channel.fromFilePairs(params.input, checkIfExists: true)
        .map { id, files -> [ [id:id], files ] }

    // 2. Quality Control & Trimming
    FASTP(ch_input)

    // 3. Alignment
    ch_fasta = file(params.fasta)
    // Collect all index constituent files into a file list channel
    ch_bwa_index = Channel
        .fromPath(params.bwa_index, checkIfExists: true)
        .collect()
    // ch_bwa_index = file(params.bwa_index)
    BWAMEM2_ALIGN(FASTP.out.reads, ch_fasta, ch_bwa_index)
    // 4. Mark Duplicates (stabilizes depth coverage stats)
    PICARD_MARKDUPLICATES(BWAMEM2_ALIGN.out.bam)



    // 5. Genotype Likelihood Generation
    // Generates un-called BCF files containing raw genotype probabilities
    BCFTOOLS_MPILEUP(PICARD_MARKDUPLICATES.out.bam_and_index, ch_fasta)

    ch_fastp_metrics  = FASTP.out.json.collect { it[1] }
    ch_picard_metrics = PICARD_MARKDUPLICATES.out.metrics.collect { it[1] }

    MULTIQC(ch_fastp_metrics, ch_picard_metrics)
    // 6. Imputation via GLIMPSE2
    // We pass the likelihoods, maps, and reference panel
    // ch_maps = file(params.genetic_map)
    // ch_ref  = file(params.ref_panel_bcf)
    // =========================================================================
    // ⚙️ INTEGRACIÓ DE LA PIPELINE WGS GLIMPSE2 (CHRS 1-22)
    // =========================================================================
    
    // Generar canal matriu per als 22 autosomes
    ch_chromosomes = Channel.of(1..22).map { "chr${it}" }

    // Preparació dels canals de referències basats en patrons parametritzats
    // Nota: params.ref_dir s'ha de passar al config o CLI
    ch_ref_sites = ch_chromosomes.map { chr -> [ chr, file("${params.ref_dir}/1000GP.${chr}.sites.vcf.gz"), file("${params.ref_dir}/1000GP.${chr}.sites.vcf.gz.csi") ] }
    ch_ref_panels = ch_chromosomes.map { chr -> [ chr, file("${params.ref_dir}/1000GP.${chr}.bcf"), file("${params.ref_dir}/1000GP.${chr}.bcf.csi") ] }
    ch_maps = ch_chromosomes.map { chr -> [ chr, file("${params.map_dir}/${chr}.b38.gmap.gz") ] }

    // PAS A: Executar generació de Chunks genòmics per cada cromosoma
    ch_chunk_input = ch_ref_sites.join(ch_maps) // Estructura: [chr, sites, tbi, map]
    GLIMPSE2_CHUNK(
        ch_chunk_input.map{ it[0] }, 
        ch_chunk_input.map{ it[1] }, 
        ch_chunk_input.map{ it[2] }, 
        ch_chunk_input.map{ it[3] }
    )

    // PAS B: Parsejar dinàmicament els fitxers de text de fragments genònics
    // Corregido: Extraemos solo el archivo del tuple [chr, file] mediante .map { it[1] }
    ch_parsed_chunks = GLIMPSE2_CHUNK.out.chunks
        .map { chr, file -> file }
        .splitText() { line ->
            def tokens = line.trim().split(/\s+/)
            return [tokens[1], tokens[0], tokens[2], tokens[3]] // [chr, chunk_id, irg, org]
        }

    // PAS C: Executar la fragmentació del panell de referència a format binari optimitzat
    // Creuem les coordenades obtingudes amb el panell i mapa genètic del seu corresponent cromosoma
    ch_split_ref_input = ch_parsed_chunks.join(ch_ref_panels).join(ch_maps)
    GLIMPSE2_SPLIT_REFERENCE(
        ch_split_ref_input.map{ [it[0], it[1], it[2], it[3]] }, // [chr, chunk_id, irg, org]
        ch_split_ref_input.map{ it[4] },                         // ref_panel bcf
        ch_split_ref_input.map{ it[5] },                         // ref_panel csi
        ch_split_ref_input.map{ it[6] }                          // genetic map
    )

    // PAS D: Llançament Massiu en Paral·lel de la Fase de Fesat / Imputació
    // Conectamos los BAM/BAI mapeados directamente con los binaris de GLIMPSE2
    ch_phase_input = PICARD_MARKDUPLICATES.out.bam_and_index // Estructura esperada: [meta, bam, bai]
        .combine(GLIMPSE2_SPLIT_REFERENCE.out.bin_ref)       // Cruce masivo: [meta, bam, bai, chr, irg, org, bin_ref]

    GLIMPSE2_PHASE(ch_phase_input)

    // PAS E: Reagrupament i Lligament (Ligate) dels fragments imputats
    ch_ligate_input = GLIMPSE2_PHASE.out.chunk_bcf
        .groupTuple(by: [0, 1]) 

    GLIMPSE2_LIGATE(ch_ligate_input)

    // El resultat final unificat de tota la genòmica del pacient es troba a:
    // GLIMPSE2_LIGATE.out.imputed_vcf
    // // 7. Polygenic Risk Score Calculation
    // // Aggregates imputed data and computes the PRS scoring matrix
    // ch_prs_weights = file(params.prs_score_file)
    // PLINK2_PRS(GLIMPSE2_PHASE.out.imputed_vcf, ch_prs_weights)
}