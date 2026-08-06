// Import raw process components
include { BCFTOOLS_MERGE_SAMPLES;
          BCFTOOLS_CONCAT_CHROME_WIDE;
          PLINK2_GWAS_QC;
          PLINK2_PCA;
          PLINK2_ASSOCIATION;
          GWAS_VISUALIZATION } from '../modules/gwas'

workflow GWAS_ANALYSIS_TRACK {
    take:
    ch_imputed_vcfs  // Channel output from GLIMPSE2_LIGATE: [meta, chr, vcf, csi]
    ch_phenotypes    // File channel pointing to phenotype text file

    main:
    // 1. Group individual files by chromosome [chr, [vcf1, vcf2...], [csi1, csi2...]]
    ch_cohort_by_chr = ch_imputed_vcfs
        .map { meta, chr, vcf, csi -> [ chr, vcf, csi ] }
        .groupTuple(by: 0)

    // 2. Perform cohort merge per chromosome
    BCFTOOLS_MERGE_SAMPLES(ch_cohort_by_chr)

    // 3. Transform channel data into structural lists crossing the chromosome boundary
    ch_vcfs_to_concat = BCFTOOLS_MERGE_SAMPLES.out.merged_vcf.map { chr, vcf, csi -> vcf }.collect()
    ch_csis_to_concat = BCFTOOLS_MERGE_SAMPLES.out.merged_vcf.map { chr, vcf, csi -> csi }.collect()

    // 4. Concat to a single genome-wide cohort VCF file
    BCFTOOLS_CONCAT_CHROME_WIDE(ch_vcfs_to_concat, ch_csis_to_concat)

    // 5. Run Quality Control & Sample Filtering
    PLINK2_GWAS_QC(BCFTOOLS_CONCAT_CHROME_WIDE.out.genome_wide_vcf)

    // 6. Run Principal Component Analysis
    PLINK2_PCA(PLINK2_GWAS_QC.out.pfile)

    // 7. Run Dosage-based Generalized Linear Model Regression
    PLINK2_ASSOCIATION(
        PLINK2_GWAS_QC.out.pfile,
        PLINK2_PCA.out.eigenvec,
        ch_phenotypes
    )

    // 8. Generate Quality Control Plots
    GWAS_VISUALIZATION(PLINK2_ASSOCIATION.out.raw_results)

    emit:
    gwas_raw_data = PLINK2_ASSOCIATION.out.raw_results
    gwas_plots    = GWAS_VISUALIZATION.out.plots
}