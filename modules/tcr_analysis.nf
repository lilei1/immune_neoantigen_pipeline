/*
========================================================================================
    TCR ANALYSIS MODULES
========================================================================================
*/

process MIXCR_ANALYZE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::mixcr=4.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mixcr:4.3.2--hdfd78af_0' :
        'mgibio/mixcr:latest' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.clns")    , emit: clns
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
        mixcr analyze shotgun \\
            -s $species \\
            --starting-material rna \\
            --only-productive \\
            --report ${prefix}.report \\
            $args \\
            $reads \\
            $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mixcr: \$(mixcr --version 2>&1 | grep -o 'MiXCR v[0-9.]*' | sed 's/MiXCR v//')
        END_VERSIONS
        """
    } else {
        """
        mixcr analyze shotgun \\
            -s $species \\
            --starting-material rna \\
            --only-productive \\
            --report ${prefix}.report \\
            $args \\
            ${reads[0]} ${reads[1]} \\
            $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mixcr: \$(mixcr --version 2>&1 | grep -o 'MiXCR v[0-9.]*' | sed 's/MiXCR v//')
        END_VERSIONS
        """
    }
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
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    mixcr exportClones \\
        --chains ALL \\
        --preset full \\
        $args \\
        $clns \\
        ${prefix}.clonotypes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixcr: \$(mixcr --version 2>&1 | grep -o 'MiXCR v[0-9.]*' | sed 's/MiXCR v//')
    END_VERSIONS
    """
}

process TRACK_CLONES {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::r-base=4.2.1 bioconda::r-ggplot2=3.3.6 conda-forge::r-dplyr=1.0.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'rocker/tidyverse:4.2.1' }"

    input:
    tuple val(meta), path(clonotype_files), val(sample_metas)

    output:
    tuple val(meta), path("*.tracking_results.tsv"), emit: tracking
    tuple val(meta), path("*.clonotype_plots.pdf") , emit: plots
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    #!/usr/bin/env Rscript
    
    library(dplyr)
    library(ggplot2)
    library(readr)
    
    # Read clonotype files
    clonotype_data <- list()
    file_names <- list.files(pattern = "*.clonotypes.txt")
    
    for (file in file_names) {
        sample_id <- gsub(".clonotypes.txt", "", file)
        data <- read_tsv(file, show_col_types = FALSE)
        data\$sample_id <- sample_id
        clonotype_data[[sample_id]] <- data
    }
    
    # Combine all data
    combined_data <- bind_rows(clonotype_data)
    
    # Track clonotypes across timepoints
    if (nrow(combined_data) > 0) {
        # Identify shared clonotypes
        clone_tracking <- combined_data %>%
            group_by(cloneId) %>%
            summarise(
                samples = paste(unique(sample_id), collapse = ";"),
                timepoints = n_distinct(sample_id),
                total_reads = sum(cloneCount, na.rm = TRUE),
                .groups = "drop"
            ) %>%
            arrange(desc(total_reads))
        
        # Save tracking results
        write_tsv(clone_tracking, "${prefix}.tracking_results.tsv")
        
        # Create plots
        pdf("${prefix}.clonotype_plots.pdf", width = 10, height = 8)
        
        # Plot 1: Clonotype diversity
        p1 <- combined_data %>%
            group_by(sample_id) %>%
            summarise(
                total_clones = n(),
                total_reads = sum(cloneCount, na.rm = TRUE),
                .groups = "drop"
            ) %>%
            ggplot(aes(x = sample_id, y = total_clones)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            theme_minimal() +
            labs(title = "Clonotype Diversity Across Samples",
                 x = "Sample", y = "Number of Clonotypes") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p1)
        
        # Plot 2: Shared clonotypes
        if (nrow(clone_tracking) > 0) {
            p2 <- clone_tracking %>%
                filter(timepoints > 1) %>%
                head(20) %>%
                ggplot(aes(x = reorder(cloneId, total_reads), y = total_reads)) +
                geom_bar(stat = "identity", fill = "coral") +
                coord_flip() +
                theme_minimal() +
                labs(title = "Top 20 Shared Clonotypes",
                     x = "Clone ID", y = "Total Read Count")
            print(p2)
        }
        
        dev.off()
    } else {
        # Create empty files if no data
        write_tsv(data.frame(), "${prefix}.tracking_results.tsv")
        pdf("${prefix}.clonotype_plots.pdf")
        plot.new()
        text(0.5, 0.5, "No clonotype data available", cex = 2)
        dev.off()
    }
    
    # Write versions
    writeLines(c(
        '"${task.process}":',
        paste0('    r-base: ', R.version.string),
        paste0('    dplyr: ', packageVersion("dplyr")),
        paste0('    ggplot2: ', packageVersion("ggplot2"))
    ), "versions.yml")
    """
}
