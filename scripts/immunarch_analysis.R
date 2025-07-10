#!/usr/bin/env Rscript
# Comprehensive TCR analysis with immunarch
# This script demonstrates diversity and clonal expansion metrics

# Install immunarch if not available
if (!require("immunarch", quietly = TRUE)) {
  if (!require("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cran.r-project.org")
  }
  devtools::install_github("immunomind/immunarch")
}

# Load required libraries
library(immunarch)
library(dplyr)
library(ggplot2)
library(readr)
library(gridExtra)
library(reshape2)

# Function to convert MiXCR output to immunarch format
convert_mixcr_to_immunarch <- function(file_path) {
  # Read MiXCR clonotype file
  data <- read_tsv(file_path, show_col_types = FALSE)
  
  # Convert to immunarch format
  # immunarch expects: CDR3.nt, CDR3.aa, V.name, J.name, Clones
  immunarch_data <- data.frame(
    CDR3.nt = data$nSeqCDR3,
    CDR3.aa = data$aaSeqCDR3,
    V.name = data$allVHitsWithScore,
    J.name = data$allJHitsWithScore,
    Clones = data$cloneCount,
    stringsAsFactors = FALSE
  )
  
  # Clean V and J gene names (take first hit)
  immunarch_data$V.name <- sapply(strsplit(immunarch_data$V.name, "\\*"), function(x) x[1])
  immunarch_data$J.name <- sapply(strsplit(immunarch_data$J.name, "\\*"), function(x) x[1])
  
  # Remove rows with missing CDR3 sequences
  immunarch_data <- immunarch_data[!is.na(immunarch_data$CDR3.aa) & 
                                   immunarch_data$CDR3.aa != "", ]
  
  return(immunarch_data)
}

# Function to analyze diversity metrics
analyze_diversity <- function(immdata) {
  cat("\n=== DIVERSITY ANALYSIS ===\n")
  
  # 1. Shannon diversity
  shannon_div <- repDiversity(immdata$data, "shannon")
  cat("Shannon diversity index:\n")
  print(shannon_div)
  
  # 2. Simpson diversity
  simpson_div <- repDiversity(immdata$data, "simpson")
  cat("\nSimpson diversity index:\n")
  print(simpson_div)
  
  # 3. Chao1 diversity
  chao1_div <- repDiversity(immdata$data, "chao1")
  cat("\nChao1 diversity index:\n")
  print(chao1_div)
  
  # 4. Hill diversity
  hill_div <- repDiversity(immdata$data, "hill")
  cat("\nHill diversity index:\n")
  print(hill_div)
  
  # 5. Rarefaction analysis
  rare <- repDiversity(immdata$data, "rarefaction")
  
  # Return all diversity metrics
  return(list(
    shannon = shannon_div,
    simpson = simpson_div,
    chao1 = chao1_div,
    hill = hill_div,
    rarefaction = rare
  ))
}

# Function to analyze clonal expansion
analyze_clonal_expansion <- function(immdata) {
  cat("\n=== CLONAL EXPANSION ANALYSIS ===\n")
  
  # 1. Clonal proportion
  clonal_props <- repClonality(immdata$data, "clonal.prop")
  cat("Clonal proportions:\n")
  print(clonal_props)
  
  # 2. Top clones
  top_clones <- repClonality(immdata$data, "top", .head = c(10, 100, 1000))
  cat("\nTop clones proportion:\n")
  print(top_clones)
  
  # 3. Homeostasis
  homeo <- repClonality(immdata$data, "homeo")
  cat("\nHomeostasis:\n")
  print(homeo)
  
  # Return all clonality metrics
  return(list(
    clonal_prop = clonal_props,
    top_clones = top_clones,
    homeostasis = homeo
  ))
}

# Function to analyze longitudinal tracking
analyze_longitudinal <- function(immdata) {
  cat("\n=== LONGITUDINAL TRACKING ===\n")
  
  # Check if we have multiple timepoints
  if (length(unique(immdata$meta$Timepoint)) <= 1) {
    cat("Not enough timepoints for longitudinal analysis\n")
    return(NULL)
  }
  
  # Track clones across timepoints
  tracking_results <- list()
  
  for (patient in unique(immdata$meta$Patient)) {
    cat(sprintf("\nTracking clones for patient: %s\n", patient))
    
    # Get samples for this patient
    patient_samples <- immdata$meta$Sample[immdata$meta$Patient == patient]
    
    if (length(patient_samples) > 1) {
      patient_data <- immdata$data[patient_samples]
      
      # Public clonotypes
      pub <- repOverlap(patient_data, "public")
      cat(sprintf("Public clonotypes: %d\n", sum(pub)))
      
      # Overlap analysis
      overlap <- repOverlap(patient_data, "overlap")
      cat("Overlap between timepoints:\n")
      print(overlap)
      
      # Track specific clones
      tracking_results[[patient]] <- list(
        public = pub,
        overlap = overlap
      )
    }
  }
  
  return(tracking_results)
}

# Function to create visualizations
create_visualizations <- function(immdata, diversity_metrics, clonality_metrics, tracking_results, output_pdf) {
  cat("\n=== CREATING VISUALIZATIONS ===\n")
  
  pdf(output_pdf, width = 12, height = 8)
  
  # 1. Diversity plots
  p1 <- vis(diversity_metrics$shannon, .by = "Patient", .meta = immdata$meta)
  print(p1 + ggtitle("Shannon Diversity Index"))
  
  p2 <- vis(diversity_metrics$simpson, .by = "Patient", .meta = immdata$meta)
  print(p2 + ggtitle("Simpson Diversity Index"))
  
  # Rarefaction curve
  p3 <- vis(diversity_metrics$rarefaction)
  print(p3 + ggtitle("Rarefaction Analysis"))
  
  # 2. Clonal expansion plots
  p4 <- vis(clonality_metrics$clonal_prop, .by = "Patient", .meta = immdata$meta)
  print(p4 + ggtitle("Clonal Proportions"))
  
  p5 <- vis(clonality_metrics$top_clones, .by = "Patient", .meta = immdata$meta)
  print(p5 + ggtitle("Top Clones Analysis"))
  
  # 3. Repertoire overlap
  if (length(immdata$data) > 1) {
    p6 <- vis(repOverlap(immdata$data, "public"))
    print(p6 + ggtitle("Public Clonotypes"))
    
    p7 <- vis(repOverlap(immdata$data, "overlap"))
    print(p7 + ggtitle("Repertoire Overlap"))
  }
  
  # 4. V/J gene usage
  gene_usage_v <- geneUsage(immdata$data, "V")
  p8 <- vis(gene_usage_v, .by = "Patient", .meta = immdata$meta)
  print(p8 + ggtitle("V Gene Usage"))
  
  gene_usage_j <- geneUsage(immdata$data, "J")
  p9 <- vis(gene_usage_j, .by = "Patient", .meta = immdata$meta)
  print(p9 + ggtitle("J Gene Usage"))
  
  # 5. Spectratype (CDR3 length distribution)
  spectra <- repExplore(immdata$data, .method = "len", .col = "aa")
  p10 <- vis(spectra)
  print(p10 + ggtitle("CDR3 Length Distribution"))
  
  # 6. Tracking plots (if available)
  if (!is.null(tracking_results)) {
    for (patient in names(tracking_results)) {
      if (!is.null(tracking_results[[patient]]$overlap)) {
        p11 <- vis(tracking_results[[patient]]$overlap)
        print(p11 + ggtitle(paste0("Clone Tracking - Patient ", patient)))
      }
    }
  }
  
  dev.off()
  cat(sprintf("Visualizations saved to: %s\n", output_pdf))
}

# Main function to run the analysis
run_immunarch_analysis <- function(input_dir, output_dir = "immunarch_results") {
  cat("=== IMMUNARCH TCR ANALYSIS ===\n")
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE)
  
  # Find all clonotype files
  clonotype_files <- list.files(input_dir, pattern = "*.clonotypes.tsv", full.names = TRUE)
  
  if (length(clonotype_files) == 0) {
    cat("No clonotype files found in:", input_dir, "\n")
    return(NULL)
  }
  
  cat("Found", length(clonotype_files), "clonotype files\n")
  
  # Process all files
  immunarch_data <- list()
  sample_names <- c()
  
  for (file in clonotype_files) {
    sample_name <- gsub("_filtered.clonotypes.tsv|_collapsed.clonotypes.tsv|.clonotypes.tsv", "", basename(file))
    sample_names <- c(sample_names, sample_name)
    
    cat("Processing:", sample_name, "\n")
    
    tryCatch({
      immunarch_data[[sample_name]] <- convert_mixcr_to_immunarch(file)
      cat("  Clonotypes:", nrow(immunarch_data[[sample_name]]), "\n")
    }, error = function(e) {
      cat("  Error processing", file, ":", e$message, "\n")
      immunarch_data[[sample_name]] <- data.frame()
    })
  }
  
  # Remove empty datasets
  immunarch_data <- immunarch_data[sapply(immunarch_data, nrow) > 0]
  
  if (length(immunarch_data) == 0) {
    cat("No valid data to analyze\n")
    return(NULL)
  }
  
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
  
  # Run analyses
  diversity_metrics <- analyze_diversity(immdata)
  clonality_metrics <- analyze_clonal_expansion(immdata)
  tracking_results <- analyze_longitudinal(immdata)
  
  # Create visualizations
  output_pdf <- file.path(output_dir, "immunarch_analysis.pdf")
  create_visualizations(immdata, diversity_metrics, clonality_metrics, tracking_results, output_pdf)
  
  # Save results
  diversity_df <- data.frame(
    Sample = names(diversity_metrics$shannon),
    Shannon = as.numeric(diversity_metrics$shannon),
    Simpson = as.numeric(diversity_metrics$simpson),
    Chao1 = as.numeric(diversity_metrics$chao1),
    Hill = as.numeric(diversity_metrics$hill),
    stringsAsFactors = FALSE
  )
  
  # Add metadata
  diversity_df <- merge(diversity_df, metadata, by = "Sample", all.x = TRUE)
  
  # Save diversity metrics
  write_tsv(diversity_df, file.path(output_dir, "diversity_metrics.tsv"))
  
  # Save clonality metrics
  clonality_df <- data.frame(
    Sample = names(clonality_metrics$clonal_prop),
    Clonal_Proportion = as.numeric(clonality_metrics$clonal_prop),
    Top10_Proportion = as.numeric(clonality_metrics$top_clones[["10"]]),
    Top100_Proportion = as.numeric(clonality_metrics$top_clones[["100"]]),
    Top1000_Proportion = as.numeric(clonality_metrics$top_clones[["1000"]]),
    stringsAsFactors = FALSE
  )
  
  # Add metadata
  clonality_df <- merge(clonality_df, metadata, by = "Sample", all.x = TRUE)
  
  # Save clonality metrics
  write_tsv(clonality_df, file.path(output_dir, "clonal_expansion.tsv"))
  
  cat("\nAnalysis complete. Results saved to:", output_dir, "\n")
  
  return(immdata)
}

# If script is run directly
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Usage: Rscript immunarch_analysis.R <input_directory> [output_directory]\n")
    quit(status = 1)
  }
  
  input_dir <- args[1]
  output_dir <- if (length(args) >= 2) args[2] else "immunarch_results"
  
  run_immunarch_analysis(input_dir, output_dir)
}
