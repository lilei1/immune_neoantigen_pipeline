#!/bin/bash

# Complete TCR Analysis Pipeline with Plots
# This script runs the full pipeline from data generation to final plots

set -e  # Exit on any error

echo "🧬 COMPLETE TCR ANALYSIS PIPELINE WITH PLOTS"
echo "=============================================="

# Step 1: Generate TCR test data (if not already done)
echo "📊 Step 1: Generating TCR longitudinal test data..."
if [ ! -d "test_data_tcr" ]; then
    python3 scripts/generate_tcr_longitudinal_data.py
    echo "✅ TCR test data generated"
else
    echo "✅ TCR test data already exists"
fi

# Step 2: Run individual sample processing with MiXCR
echo "📊 Step 2: Processing samples with MiXCR..."

# Create output directory
mkdir -p tcr_analysis_results

# Process each sample
samples=(
    "PATIENT_01_T0_baseline"
    "PATIENT_01_T1_cycle1" 
    "PATIENT_01_T2_cycle3"
    "PATIENT_01_T3_progression"
    "PATIENT_01_T4_posttreatment"
    "PATIENT_02_T0_baseline"
    "PATIENT_02_T1_cycle1"
    "PATIENT_02_T2_cycle3" 
    "PATIENT_02_T3_progression"
    "PATIENT_02_T4_posttreatment"
)

for sample in "${samples[@]}"; do
    echo "  Processing $sample..."
    
    # Check if files exist
    r1_file="test_data_tcr/${sample}_tcr_1.fastq.gz"
    r2_file="test_data_tcr/${sample}_tcr_2.fastq.gz"
    
    if [ ! -f "$r1_file" ] || [ ! -f "$r2_file" ]; then
        echo "    ⚠️  Skipping $sample - files not found"
        continue
    fi
    
    # MiXCR processing
    output_dir="tcr_analysis_results/$sample"
    mkdir -p "$output_dir"
    
    # Align
    if [ ! -f "$output_dir/${sample}.vdjca" ]; then
        echo "    🔄 Aligning..."
        docker run --rm -v $(pwd):/data -w /data mgibio/mixcr:latest \
            mixcr align -s hsa \
            --report "/data/$output_dir/${sample}_align.report" \
            "/data/$r1_file" "/data/$r2_file" \
            "/data/$output_dir/${sample}.vdjca"
    fi
    
    # Assemble
    if [ ! -f "$output_dir/${sample}.clns" ]; then
        echo "    🔄 Assembling..."
        docker run --rm -v $(pwd):/data -w /data mgibio/mixcr:latest \
            mixcr assemble \
            --report "/data/$output_dir/${sample}_assemble.report" \
            "/data/$output_dir/${sample}.vdjca" \
            "/data/$output_dir/${sample}.clns"
    fi
    
    # Export clonotypes
    if [ ! -f "$output_dir/${sample}.clonotypes.tsv" ]; then
        echo "    🔄 Exporting clonotypes..."
        docker run --rm -v $(pwd):/data -w /data mgibio/mixcr:latest \
            mixcr exportClones --chains ALL --preset full -f \
            "/data/$output_dir/${sample}.clns" \
            "/data/$output_dir/${sample}.clonotypes.tsv"
    fi
    
    echo "    ✅ $sample completed"
done

echo "✅ MiXCR processing completed"

# Step 3: Run immunarch analysis to generate plots
echo "📊 Step 3: Running immunarch analysis and generating plots..."

# Create a directory with all clonotype files for immunarch
mkdir -p tcr_analysis_results/all_clonotypes
for sample in "${samples[@]}"; do
    clonotype_file="tcr_analysis_results/$sample/${sample}.clonotypes.tsv"
    if [ -f "$clonotype_file" ]; then
        cp "$clonotype_file" "tcr_analysis_results/all_clonotypes/"
    fi
done

# Run immunarch analysis
Rscript scripts/immunarch_analysis.R tcr_analysis_results/all_clonotypes tcr_analysis_results/immunarch_output

echo "✅ Immunarch analysis completed"

# Step 4: Create summary plots with Python
echo "📊 Step 4: Creating additional summary plots..."

python3 - << 'EOF'
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import glob

# Set style
plt.style.use('default')
sns.set_palette("husl")

# Read diversity metrics
diversity_file = "tcr_analysis_results/immunarch_output/diversity_metrics.tsv"
expansion_file = "tcr_analysis_results/immunarch_output/clonal_expansion.tsv"

output_dir = Path("tcr_analysis_results/summary_plots")
output_dir.mkdir(exist_ok=True)

try:
    # Load data
    diversity_df = pd.read_csv(diversity_file, sep='\t')
    expansion_df = pd.read_csv(expansion_file, sep='\t')
    
    # Create summary plots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('TCR Repertoire Analysis - 5 Timepoint Longitudinal Study', fontsize=16, fontweight='bold')
    
    # Plot 1: Shannon Diversity by Patient and Timepoint
    ax1 = axes[0, 0]
    for patient in diversity_df['Patient'].unique():
        patient_data = diversity_df[diversity_df['Patient'] == patient]
        ax1.plot(range(len(patient_data)), patient_data['Shannon'], 'o-', label=patient, linewidth=2, markersize=8)
    ax1.set_title('Shannon Diversity Over Time', fontweight='bold')
    ax1.set_xlabel('Timepoint')
    ax1.set_ylabel('Shannon Diversity Index')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Simpson Diversity
    ax2 = axes[0, 1]
    for patient in diversity_df['Patient'].unique():
        patient_data = diversity_df[diversity_df['Patient'] == patient]
        ax2.plot(range(len(patient_data)), patient_data['Simpson'], 's-', label=patient, linewidth=2, markersize=8)
    ax2.set_title('Simpson Diversity Over Time', fontweight='bold')
    ax2.set_xlabel('Timepoint')
    ax2.set_ylabel('Simpson Diversity Index')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Chao1 Diversity
    ax3 = axes[0, 2]
    for patient in diversity_df['Patient'].unique():
        patient_data = diversity_df[diversity_df['Patient'] == patient]
        ax3.plot(range(len(patient_data)), patient_data['Chao1'], '^-', label=patient, linewidth=2, markersize=8)
    ax3.set_title('Chao1 Diversity Over Time', fontweight='bold')
    ax3.set_xlabel('Timepoint')
    ax3.set_ylabel('Chao1 Diversity Index')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Clonal Proportion
    ax4 = axes[1, 0]
    for patient in expansion_df['Patient'].unique():
        patient_data = expansion_df[expansion_df['Patient'] == patient]
        ax4.plot(range(len(patient_data)), patient_data['Clonal_Proportion'], 'o-', label=patient, linewidth=2, markersize=8)
    ax4.set_title('Clonal Proportion Over Time', fontweight='bold')
    ax4.set_xlabel('Timepoint')
    ax4.set_ylabel('Clonal Proportion')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: Top 10 Clone Proportion
    ax5 = axes[1, 1]
    for patient in expansion_df['Patient'].unique():
        patient_data = expansion_df[expansion_df['Patient'] == patient]
        ax5.plot(range(len(patient_data)), patient_data['Top10_Proportion'], 's-', label=patient, linewidth=2, markersize=8)
    ax5.set_title('Top 10 Clones Proportion Over Time', fontweight='bold')
    ax5.set_xlabel('Timepoint')
    ax5.set_ylabel('Top 10 Clones Proportion')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: Diversity vs Clonal Expansion
    ax6 = axes[1, 2]
    scatter = ax6.scatter(diversity_df['Shannon'], expansion_df['Clonal_Proportion'], 
                         c=range(len(diversity_df)), cmap='viridis', s=100, alpha=0.7)
    ax6.set_title('Diversity vs Clonal Expansion', fontweight='bold')
    ax6.set_xlabel('Shannon Diversity')
    ax6.set_ylabel('Clonal Proportion')
    plt.colorbar(scatter, ax=ax6, label='Sample Order')
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'tcr_longitudinal_summary.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'tcr_longitudinal_summary.pdf', bbox_inches='tight')
    
    print("✅ Summary plots created successfully!")
    print(f"📁 Plots saved to: {output_dir}")
    
except Exception as e:
    print(f"❌ Error creating plots: {e}")
    print("This might be because immunarch analysis hasn't completed yet")

EOF

echo "✅ Additional summary plots completed"

# Step 5: Display results summary
echo ""
echo "🎉 PIPELINE COMPLETED SUCCESSFULLY!"
echo "=================================="
echo ""
echo "📁 Results are available in:"
echo "   • tcr_analysis_results/immunarch_output/immunarch_analysis.pdf"
echo "   • tcr_analysis_results/summary_plots/tcr_longitudinal_summary.png"
echo "   • tcr_analysis_results/summary_plots/tcr_longitudinal_summary.pdf"
echo ""
echo "📊 Data files:"
echo "   • tcr_analysis_results/immunarch_output/diversity_metrics.tsv"
echo "   • tcr_analysis_results/immunarch_output/clonal_expansion.tsv"
echo ""
echo "🔍 Individual sample results:"
for sample in "${samples[@]}"; do
    if [ -d "tcr_analysis_results/$sample" ]; then
        echo "   • tcr_analysis_results/$sample/"
    fi
done
echo ""
echo "🚀 To view plots:"
echo "   open tcr_analysis_results/summary_plots/tcr_longitudinal_summary.pdf"
echo "   open tcr_analysis_results/immunarch_output/immunarch_analysis.pdf"
