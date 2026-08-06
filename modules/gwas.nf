process BCFTOOLS_MERGE_SAMPLES {
    tag "${chr}"
    label 'process_high'

    input:
    tuple val(chr), path(vcfs), path(csis)

    output:
    tuple val(chr), path("${chr}_cohort_merged.vcf.gz"), path("${chr}_cohort_merged.vcf.gz.csi"), emit: merged_chr_vcf

    script:
    """
    bcftools merge \\
        --threads ${task.cpus} \\
        -O z \\
        -o ${chr}_cohort_merged.vcf.gz \\
        ${vcfs}

    bcftools index \\
        --threads ${task.cpus} \\
        ${chr}_cohort_merged.vcf.gz
    """
}

process BCFTOOLS_CONCAT_CHROME_WIDE {
    label 'process_high'

    input:
    path(vcfs)
    path(csis)

    output:
    tuple path("cohort_genome_wide.vcf.gz"), path("cohort_genome_wide.vcf.gz.csi"), emit: genome_wide_vcf

    script:
    // Ensure chromosomes are sorted numerically to preserve coordinate integrity
    def sorted_vcfs = vcfs.sort { it.name }
    """
    bcftools concat \\
        --threads ${task.cpus} \\
        -O z \\
        -o cohort_genome_wide.vcf.gz \\
        ${sorted_vcfs}

    bcftools index \\
        --threads ${task.cpus} \\
        cohort_genome_wide.vcf.gz
    """
}

process PLINK2_GWAS_QC {
    label 'process_high'

    input:
    tuple path(vcf), path(csi)

    output:
    tuple path("cohort_clean.pgen"), path("cohort_clean.pvar"), path("cohort_clean.psam"), emit: pfile

    script:
    """
    plink2 --vcf ${vcf} \\
           --vcf-dosage DS \\
           --maf 0.01 \\
           --hwe 1e-6 \\
           --geno 0.02 \\
           --make-pgen \\
           --out cohort_clean
    """
}

process PLINK2_PCA {
    label 'process_high'

    input:
    tuple path(pgen), path(pvar), path(psam)

    output:
    path("cohort_pca.eigenvec"), emit: eigenvec
    path("cohort_pca.eigenval"), emit: eigenval

    script:
    """
    # 1. Variant Pruning based on Linkage Disequilibrium (LD)
    plink2 --pfile cohort_clean \\
           --indep-pairwise 50 5 0.2 \\
           --out pruned_variants

    # 2. Extract independent structural variants and run PCA
    plink2 --pfile cohort_clean \\
           --extract pruned_variants.prune.in \\
           --pca 10 \\
           --out cohort_pca
    """
}

process PLINK2_ASSOCIATION {
    label 'process_high'

    input:
    tuple path(pgen), path(pvar), path(psam)
    path eigenvec
    path phenotypes

    output:
    path("gwas_results.*"), emit: raw_results

    script:
    """
    plink2 --pfile cohort_clean \\
           --vcf-dosage DS \\
           --glm cols=chrom,pos,ref,alt,ax,a1freq,nobs,beta,orbeta,se,tz,p allow-no-covars \\
           --pheno ${phenotypes} \\
           --covar ${eigenvec} \\
           --covar-name PC1-PC4 \\
           --out gwas_results
    """
}

process GWAS_VISUALIZATION {
    label 'process_low'

    input:
    path gwas_files

    output:
    path("gwas_manhattan_qq_plots.png"), emit: plots

    script:
    """
    #!/usr/bin/env Rscript
    gwas_file <- list.files(pattern = "gwas_results\\\\..*\\\\.glm\\\\.(linear|logistic)")
    if (length(gwas_file) == 0) stop("No association outcome file detected!")
    
    data <- read.table(gwas_file[1], header=TRUE, sep="\\t", comment.char="", check.names=FALSE)
    data <- data[data\$TEST == "ADD", ]
    
    df <- data[, c("#CHROM", "POS", "P")]
    colnames(df) <- c("CHR", "BP", "P")
    df\$CHR <- as.numeric(gsub("chr", "", df\$CHR))
    df <- df[!is.na(df\$P) & !is.na(df\$CHR) & !is.na(df\$BP), ]
    
    png("gwas_manhattan_qq_plots.png", width=1400, height=600, res=150)
    par(mfrow=c(1,2))
    
    # Manhattan Plot
    df\$pos_cum <- df\$BP + (df\$CHR * 3e8)
    colors <- rep(c("#2c3e50", "#16a085"), 11)
    plot(df\$pos_cum, -log10(df\$P), col=colors[df\$CHR], pch=20, cex=0.6,
         xlab="Chromosomes (1-22)", ylab="-log10(P-value)", main="Genome-Wide Association Manhattan Plot", xaxt="n")
    abline(h=-log10(5e-8), col="#e74c3c", lty=2, lwd=1.5)
    
    # Q-Q Plot
    observed <- -log10(sort(df\$P))
    expected <- -log10(ppoints(length(df\$P)))
    plot(expected, observed, pch=20, cex=0.6, col="#2980b9",
         xlab="Expected -log10(P)", ylab="Observed -log10(P)", main="GWAS Q-Q Diagnostics Plot")
    abline(0, 1, col="#e74c3c", lty=2, lwd=1.5)
    
    dev.off()
    """
}