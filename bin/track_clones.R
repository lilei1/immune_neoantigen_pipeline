#!/usr/bin/env Rscript

# Enhanced script for clonotype tracking and visualization
# Author: lilei
# Description: Tracks T-cell clonotypes across longitudinal samples

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(tidyr)
    library(stringr)
    library(RColorBrewer)
})

# Function to read and process clonotype data
read_clonotype_data <- function(file_path, sample_id) {
    tryCatch({
        data <- read_tsv(file_path, show_col_types = FALSE)
        data$sample_id <- sample_id
        return(data)
    }, error = function(e) {
        warning(paste("Could not read file:", file_path, "Error:", e$message))
        return(NULL)
    })
}

# Function to track clonotypes across timepoints
track_clonotypes <- function(clonotype_data) {
    if (nrow(clonotype_data) == 0) {
        return(data.frame())
    }

    # Identify shared clonotypes
    clone_tracking <- clonotype_data %>%
        group_by(cloneId) %>%
        summarise(
            samples = paste(unique(sample_id), collapse = ";"),
            timepoints = n_distinct(sample_id),
            total_reads = sum(cloneCount, na.rm = TRUE),
            avg_frequency = mean(cloneFraction, na.rm = TRUE),
            max_frequency = max(cloneFraction, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        arrange(desc(total_reads))

    return(clone_tracking)
}

# Function to create diversity plots
plot_diversity <- function(clonotype_data) {
    diversity_stats <- clonotype_data %>%
        group_by(sample_id) %>%
        summarise(
            total_clones = n(),
            total_reads = sum(cloneCount, na.rm = TRUE),
            shannon_diversity = -sum(cloneFraction * log(cloneFraction), na.rm = TRUE),
            simpson_diversity = 1 - sum(cloneFraction^2, na.rm = TRUE),
            .groups = "drop"
        )

    p1 <- ggplot(diversity_stats, aes(x = sample_id, y = total_clones)) +
        geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
        theme_minimal() +
        labs(title = "Clonotype Diversity Across Samples",
             x = "Sample", y = "Number of Clonotypes") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    p2 <- ggplot(diversity_stats, aes(x = sample_id, y = shannon_diversity)) +
        geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
        theme_minimal() +
        labs(title = "Shannon Diversity Index",
             x = "Sample", y = "Shannon Diversity") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    return(list(clonotype_count = p1, shannon_diversity = p2))
}

# Function to create tracking plots
plot_tracking <- function(clone_tracking, clonotype_data) {
    if (nrow(clone_tracking) == 0) {
        return(list())
    }

    # Plot shared clonotypes
    shared_clones <- clone_tracking %>%
        filter(timepoints > 1) %>%
        head(20)

    p1 <- ggplot(shared_clones, aes(x = reorder(cloneId, total_reads), y = total_reads)) +
        geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
        coord_flip() +
        theme_minimal() +
        labs(title = "Top 20 Shared Clonotypes",
             x = "Clone ID", y = "Total Read Count")

    # Clonotype persistence plot
    persistence_data <- clone_tracking %>%
        count(timepoints, name = "clone_count") %>%
        mutate(persistence_category = case_when(
            timepoints == 1 ~ "Unique",
            timepoints == 2 ~ "Shared (2 timepoints)",
            timepoints >= 3 ~ "Persistent (3+ timepoints)"
        ))

    p2 <- ggplot(persistence_data, aes(x = persistence_category, y = clone_count)) +
        geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
        theme_minimal() +
        labs(title = "Clonotype Persistence",
             x = "Persistence Category", y = "Number of Clonotypes") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    return(list(shared_clones = p1, persistence = p2))
}

# Main execution
main <- function() {
    # Get command line arguments
    args <- commandArgs(trailingOnly = TRUE)

    if (length(args) < 2) {
        cat("Usage: track_clones.R <output_prefix> <clonotype_file1> [clonotype_file2] ...\n")
        quit(status = 1)
    }

    output_prefix <- args[1]
    clonotype_files <- args[-1]

    # Read all clonotype files
    all_data <- list()
    for (file in clonotype_files) {
        sample_id <- gsub(".clonotypes.txt", "", basename(file))
        data <- read_clonotype_data(file, sample_id)
        if (!is.null(data)) {
            all_data[[sample_id]] <- data
        }
    }

    if (length(all_data) == 0) {
        cat("No valid clonotype data found\n")
        quit(status = 1)
    }

    # Combine all data
    combined_data <- bind_rows(all_data)

    # Track clonotypes
    clone_tracking <- track_clonotypes(combined_data)

    # Save tracking results
    write_tsv(clone_tracking, paste0(output_prefix, ".tracking_results.tsv"))

    # Create plots
    pdf(paste0(output_prefix, ".clonotype_plots.pdf"), width = 12, height = 8)

    # Diversity plots
    diversity_plots <- plot_diversity(combined_data)
    print(diversity_plots$clonotype_count)
    print(diversity_plots$shannon_diversity)

    # Tracking plots
    tracking_plots <- plot_tracking(clone_tracking, combined_data)
    if (length(tracking_plots) > 0) {
        print(tracking_plots$shared_clones)
        print(tracking_plots$persistence)
    }

    dev.off()

    cat("Clonotype tracking completed successfully\n")
    cat("Results saved to:", paste0(output_prefix, ".tracking_results.tsv"), "\n")
    cat("Plots saved to:", paste0(output_prefix, ".clonotype_plots.pdf"), "\n")
}

# Run main function
if (!interactive()) {
    main()
}