/*
========================================================================================
    TCR ANALYSIS MODULES
========================================================================================
*/

process MIXCR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::mixcr=4.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mixcr:4.3.2--hdfd78af_0' :
        'mgibio/mixcr:latest' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.vdjca")   , emit: vdjca
    tuple val(meta), path("*.report")  , emit: report
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def species = params.mixcr_species ?: 'hsa'

    if (meta.single_end) {
        """
        mixcr align \\
            -s $species \\
            --report ${prefix}_align.report \\
            $args \\
            $reads \\
            ${prefix}.vdjca

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mixcr: \$(mixcr --version 2>&1 | grep -o 'MiXCR v[0-9.]*' | sed 's/MiXCR v//')
        END_VERSIONS
        """
    } else {
        """
        mixcr align \\
            -s $species \\
            --report ${prefix}_align.report \\
            $args \\
            ${reads[0]} ${reads[1]} \\
            ${prefix}.vdjca

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mixcr: \$(mixcr --version 2>&1 | grep -o 'MiXCR v[0-9.]*' | sed 's/MiXCR v//')
        END_VERSIONS
        """
    }
}

process MIXCR_ASSEMBLE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::mixcr=4.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mixcr:4.3.2--hdfd78af_0' :
        'mgibio/mixcr:latest' }"

    input:
    tuple val(meta), path(vdjca)

    output:
    tuple val(meta), path("*.clns")    , emit: clns
    tuple val(meta), path("*.report")  , emit: report
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mixcr assemble \\
        --report ${prefix}_assemble.report \\
        $args \\
        $vdjca \\
        ${prefix}.clns

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixcr: \$(mixcr --version 2>&1 | grep -o 'MiXCR v[0-9.]*' | sed 's/MiXCR v//')
    END_VERSIONS
    """
}

process MIXCR_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mixcr=4.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mixcr:4.3.2--hdfd78af_0' :
        'mgibio/mixcr:latest' }"

    input:
    tuple val(meta), path(clns)

    output:
    tuple val(meta), path("*_filtered.clns"), emit: clns_filtered
    tuple val(meta), path("*.report")       , emit: report
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_count = params.mixcr_min_clone_count ?: 2
    def min_fraction = params.mixcr_min_clone_fraction ?: 0.00001

    """
    # Filter clones by count and frequency
    mixcr filter \\
        --min-count $min_count \\
        --min-fraction $min_fraction \\
        --report ${prefix}_filter.report \\
        $args \\
        $clns \\
        ${prefix}_filtered.clns

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixcr: \$(mixcr --version 2>&1 | grep -o 'MiXCR v[0-9.]*' | sed 's/MiXCR v//')
    END_VERSIONS
    """
}

process MIXCR_COLLAPSE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::mixcr=4.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mixcr:4.3.2--hdfd78af_0' :
        'mgibio/mixcr:latest' }"

    input:
    tuple val(meta), path(clns)

    output:
    tuple val(meta), path("*_collapsed.clns"), emit: clns_collapsed
    tuple val(meta), path("*.report")        , emit: report
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def collapse_by = params.mixcr_collapse_by ?: 'CDR3'
    def collapse_threshold = params.mixcr_collapse_threshold ?: 0.9

    """
    # Collapse similar clonotypes to account for sequencing errors and somatic hypermutation
    mixcr collapse \\
        --by-$collapse_by \\
        --level $collapse_threshold \\
        --report ${prefix}_collapse.report \\
        $args \\
        $clns \\
        ${prefix}_collapsed.clns

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixcr: \$(mixcr --version 2>&1 | grep -o 'MiXCR v[0-9.]*' | sed 's/MiXCR v//')
    END_VERSIONS
    """
}

process MIXCR_EXPORTCLONES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mixcr=4.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mixcr:4.3.2--hdfd78af_0' :
        'mgibio/mixcr:latest' }"

    input:
    tuple val(meta), path(clns)

    output:
    tuple val(meta), path("*.clonotypes.txt"), emit: clonotypes
    tuple val(meta), path("*.clonotypes.tsv"), emit: clonotypes_tsv
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Export clones in multiple formats for downstream analysis
    mixcr exportClones \\
        --chains ALL \\
        --preset full \\
        $args \\
        $clns \\
        ${prefix}.clonotypes.txt

    # Export in TSV format for R analysis
    mixcr exportClones \\
        --chains ALL \\
        --preset full \\
        -f \\
        $args \\
        $clns \\
        ${prefix}.clonotypes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixcr: \$(mixcr --version 2>&1 | grep -o 'MiXCR v[0-9.]*' | sed 's/MiXCR v//')
    END_VERSIONS
    """
}

process IMMUNARCH_ANALYSIS {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "conda-forge::r-base=4.2.1 conda-forge::r-devtools=2.4.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'rocker/tidyverse:4.2.1' }"

    input:
    tuple val(meta), path(clonotype_files)

    output:
    tuple val(meta), path("*_diversity_metrics.tsv")    , emit: diversity
    tuple val(meta), path("*_clonal_expansion.tsv")     , emit: expansion
    tuple val(meta), path("*_longitudinal_tracking.tsv"), emit: tracking
    tuple val(meta), path("*_immunarch_plots.pdf")      , emit: plots
    tuple val(meta), path("*_immunarch_report.html")    , emit: report
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    # Install immunarch if not available
    if (!require("immunarch", quietly = TRUE)) {
        if (!require("devtools", quietly = TRUE)) {
            install.packages("devtools", repos = "https://cran.r-project.org")
        }
        devtools::install_github("immunomind/immunarch")
    }

    library(immunarch)
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(rmarkdown)

    # Function to convert MiXCR output to immunarch format
    convert_mixcr_to_immunarch <- function(file_path) {
        # Read MiXCR clonotype file
        data <- read_tsv(file_path, show_col_types = FALSE)

        # Convert to immunarch format
        # immunarch expects: CDR3.nt, CDR3.aa, V.name, J.name, Clones
        immunarch_data <- data.frame(
            CDR3.nt = data\$nSeqCDR3,
            CDR3.aa = data\$aaSeqCDR3,
            V.name = data\$allVHitsWithScore,
            J.name = data\$allJHitsWithScore,
            Clones = data\$cloneCount,
            stringsAsFactors = FALSE
        )

        # Clean V and J gene names (take first hit)
        immunarch_data\$V.name <- sapply(strsplit(immunarch_data\$V.name, "\\\\*"), function(x) x[1])
        immunarch_data\$J.name <- sapply(strsplit(immunarch_data\$J.name, "\\\\*"), function(x) x[1])

        # Remove rows with missing CDR3 sequences
        immunarch_data <- immunarch_data[!is.na(immunarch_data\$CDR3.aa) &
                                       immunarch_data\$CDR3.aa != "", ]

        return(immunarch_data)
    }

    # Process all clonotype files
    clonotype_files <- list.files(pattern = "*.clonotypes.tsv", full.names = TRUE)

    if (length(clonotype_files) == 0) {
        # Create empty output files
        write_tsv(data.frame(), "${prefix}_diversity_metrics.tsv")
        write_tsv(data.frame(), "${prefix}_clonal_expansion.tsv")
        write_tsv(data.frame(), "${prefix}_longitudinal_tracking.tsv")

        pdf("${prefix}_immunarch_plots.pdf")
        plot.new()
        text(0.5, 0.5, "No clonotype data available", cex = 2)
        dev.off()

        writeLines("No data available for analysis", "${prefix}_immunarch_report.html")

    } else {
        # Convert files to immunarch format
        immunarch_data <- list()
        sample_names <- c()

        for (file in clonotype_files) {
            sample_name <- gsub("_filtered.clonotypes.tsv|.clonotypes.tsv", "", basename(file))
            sample_names <- c(sample_names, sample_name)

            tryCatch({
                immunarch_data[[sample_name]] <- convert_mixcr_to_immunarch(file)
            }, error = function(e) {
                cat("Error processing", file, ":", e\$message, "\\n")
                immunarch_data[[sample_name]] <- data.frame()
            })
        }

        # Remove empty datasets
        immunarch_data <- immunarch_data[sapply(immunarch_data, nrow) > 0]

        if (length(immunarch_data) > 0) {
            # Create metadata
            metadata <- data.frame(
                Sample = names(immunarch_data),
                Patient = sapply(strsplit(names(immunarch_data), "_"), function(x) x[1]),
                Timepoint = sapply(strsplit(names(immunarch_data), "_"), function(x)
                    if(length(x) > 1) paste(x[2:length(x)], collapse = "_") else "Unknown"),
                stringsAsFactors = FALSE
            )

            # Create immunarch object
            immdata <- list(data = immunarch_data, meta = metadata)

            # 1. DIVERSITY ANALYSIS
            cat("Computing diversity metrics...\\n")

            # Shannon diversity
            shannon_div <- repDiversity(immdata\$data, "shannon")

            # Simpson diversity
            simpson_div <- repDiversity(immdata\$data, "simpson")

            # Chao1 diversity
            chao1_div <- repDiversity(immdata\$data, "chao1")

            # Hill diversity
            hill_div <- repDiversity(immdata\$data, "hill")

            # Combine diversity metrics
            diversity_metrics <- data.frame(
                Sample = names(shannon_div),
                Shannon = as.numeric(shannon_div),
                Simpson = as.numeric(simpson_div),
                Chao1 = as.numeric(chao1_div),
                Hill = as.numeric(hill_div),
                stringsAsFactors = FALSE
            )

            # Add metadata
            diversity_metrics <- merge(diversity_metrics, metadata, by = "Sample", all.x = TRUE)

            write_tsv(diversity_metrics, "${prefix}_diversity_metrics.tsv")

            # 2. CLONAL EXPANSION ANALYSIS
            cat("Analyzing clonal expansion...\\n")

            # Clonal proportion analysis
            clonal_props <- repClonality(immdata\$data, "clonal.prop")

            # Top clones analysis
            top_clones <- repClonality(immdata\$data, "top", .head = c(10, 100, 1000))

            # Combine expansion metrics
            expansion_metrics <- data.frame(
                Sample = names(clonal_props),
                Clonal_Proportion = as.numeric(clonal_props),
                Top10_Proportion = as.numeric(top_clones[["10"]]),
                Top100_Proportion = as.numeric(top_clones[["100"]]),
                Top1000_Proportion = as.numeric(top_clones[["1000"]]),
                stringsAsFactors = FALSE
            )

            # Add metadata
            expansion_metrics <- merge(expansion_metrics, metadata, by = "Sample", all.x = TRUE)

            write_tsv(expansion_metrics, "${prefix}_clonal_expansion.tsv")

            # 3. LONGITUDINAL TRACKING
            cat("Performing longitudinal tracking...\\n")

            if (length(unique(metadata\$Patient)) > 1 || length(unique(metadata\$Timepoint)) > 1) {
                # Track clones across timepoints
                tracking_results <- list()

                for (patient in unique(metadata\$Patient)) {
                    patient_samples <- metadata[metadata\$Patient == patient, "Sample"]

                    if (length(patient_samples) > 1) {
                        patient_data <- immdata\$data[patient_samples]

                        # Find shared clones
                        all_clones <- unique(unlist(lapply(patient_data, function(x) x\$CDR3.aa)))

                        clone_tracking <- data.frame(
                            Patient = patient,
                            CDR3_aa = all_clones,
                            stringsAsFactors = FALSE
                        )

                        # Add presence/absence and counts for each sample
                        for (sample in patient_samples) {
                            sample_data <- patient_data[[sample]]
                            clone_tracking[[paste0(sample, "_present")]] <-
                                clone_tracking\$CDR3_aa %in% sample_data\$CDR3.aa

                            clone_tracking[[paste0(sample, "_count")]] <-
                                sapply(clone_tracking\$CDR3_aa, function(cdr3) {
                                    idx <- which(sample_data\$CDR3.aa == cdr3)
                                    if (length(idx) > 0) sample_data\$Clones[idx[1]] else 0
                                })
                        }

                        tracking_results[[patient]] <- clone_tracking
                    }
                }

                if (length(tracking_results) > 0) {
                    combined_tracking <- do.call(rbind, tracking_results)
                    write_tsv(combined_tracking, "${prefix}_longitudinal_tracking.tsv")
                } else {
                    write_tsv(data.frame(), "${prefix}_longitudinal_tracking.tsv")
                }
            } else {
                write_tsv(data.frame(), "${prefix}_longitudinal_tracking.tsv")
            }

            # 4. VISUALIZATION
            cat("Creating visualizations...\\n")

            pdf("${prefix}_immunarch_plots.pdf", width = 12, height = 8)

            # Plot 1: Diversity comparison
            if (nrow(diversity_metrics) > 1) {
                p1 <- vis(shannon_div, .by = "Patient", .meta = metadata)
                print(p1)

                p2 <- vis(simpson_div, .by = "Patient", .meta = metadata)
                print(p2)
            }

            # Plot 2: Clonal expansion
            if (nrow(expansion_metrics) > 1) {
                p3 <- vis(clonal_props, .by = "Patient", .meta = metadata)
                print(p3)

                p4 <- vis(top_clones, .by = "Patient", .meta = metadata)
                print(p4)
            }

            # Plot 3: Repertoire overlap (if multiple samples)
            if (length(immdata\$data) > 1) {
                overlap <- repOverlap(immdata\$data, "public", .verbose = FALSE)
                p5 <- vis(overlap)
                print(p5)
            }

            # Plot 4: Gene usage
            gene_usage <- geneUsage(immdata\$data, "V")
            p6 <- vis(gene_usage, .by = "Patient", .meta = metadata)
            print(p6)

            dev.off()

            # 5. GENERATE HTML REPORT
            cat("Generating HTML report...\\n")

            report_content <- paste0(
                "<html><head><title>TCR Immunarch Analysis Report</title></head><body>",
                "<h1>TCR Repertoire Analysis Report</h1>",
                "<h2>Summary</h2>",
                "<p>Total samples analyzed: ", length(immdata\$data), "</p>",
                "<p>Total patients: ", length(unique(metadata\$Patient)), "</p>",
                "<p>Total timepoints: ", length(unique(metadata\$Timepoint)), "</p>",
                "<h2>Diversity Metrics</h2>",
                "<p>Mean Shannon diversity: ", round(mean(diversity_metrics\$Shannon, na.rm = TRUE), 3), "</p>",
                "<p>Mean Simpson diversity: ", round(mean(diversity_metrics\$Simpson, na.rm = TRUE), 3), "</p>",
                "<h2>Clonal Expansion</h2>",
                "<p>Mean clonal proportion: ", round(mean(expansion_metrics\$Clonal_Proportion, na.rm = TRUE), 3), "</p>",
                "<p>Mean top 10 clone proportion: ", round(mean(expansion_metrics\$Top10_Proportion, na.rm = TRUE), 3), "</p>",
                "<p>Analysis completed successfully. See accompanying TSV files and PDF plots for detailed results.</p>",
                "</body></html>"
            )

            writeLines(report_content, "${prefix}_immunarch_report.html")

        } else {
            # No valid data
            write_tsv(data.frame(), "${prefix}_diversity_metrics.tsv")
            write_tsv(data.frame(), "${prefix}_clonal_expansion.tsv")
            write_tsv(data.frame(), "${prefix}_longitudinal_tracking.tsv")

            pdf("${prefix}_immunarch_plots.pdf")
            plot.new()
            text(0.5, 0.5, "No valid clonotype data for analysis", cex = 2)
            dev.off()

            writeLines("No valid data available for analysis", "${prefix}_immunarch_report.html")
        }
    }

    # Write versions
    writeLines(c(
        '"${task.process}":',
        paste0('    r-base: ', R.version.string),
        '    immunarch: latest',
        paste0('    dplyr: ', packageVersion("dplyr")),
        paste0('    ggplot2: ', packageVersion("ggplot2"))
    ), "versions.yml")
    """
}
