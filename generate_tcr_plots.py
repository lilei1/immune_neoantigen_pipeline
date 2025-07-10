#!/usr/bin/env python3
"""
Generate TCR analysis plots directly from MiXCR output files.
This script creates diversity and clonal expansion plots for longitudinal TCR data.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import glob
import re
import random

def extract_sample_info(filename):
    """Extract patient and timepoint from filename."""
    basename = os.path.basename(filename)
    match = re.match(r'(PATIENT_\d+)_(T\d+_\w+)', basename)
    if match:
        return match.group(1), match.group(2)
    else:
        parts = basename.split('_')
        if len(parts) >= 2:
            return parts[0], '_'.join(parts[1:])
        return basename, "unknown"

def load_clonotype_data(directory):
    """Load all clonotype files from a directory."""
    clonotype_files = glob.glob(os.path.join(directory, "*/*.clonotypes.tsv"))
    
    if not clonotype_files:
        print(f"No clonotype files found in {directory}")
        return None
    
    print(f"Found {len(clonotype_files)} clonotype files")
    
    # Process each file
    all_data = []
    for file_path in clonotype_files:
        try:
            # Read the file
            df = pd.read_csv(file_path, sep='\t')
            
            # Extract sample info
            patient, timepoint = extract_sample_info(file_path)
            
            # Add metadata
            df['patient'] = patient
            df['timepoint'] = timepoint
            df['sample'] = f"{patient}_{timepoint}"
            df['file'] = os.path.basename(file_path)
            
            # Calculate basic metrics
            total_clones = len(df)
            total_reads = df['cloneCount'].sum() if 'cloneCount' in df.columns else 0
            
            # Create a summary row
            summary = {
                'patient': patient,
                'timepoint': timepoint,
                'sample': f"{patient}_{timepoint}",
                'total_clones': total_clones,
                'total_reads': total_reads,
                'file': os.path.basename(file_path)
            }
            
            # Add diversity metrics
            if 'cloneCount' in df.columns and total_clones > 0:
                # Shannon diversity
                if total_reads > 0:
                    frequencies = df['cloneCount'] / total_reads
                    shannon = -sum(frequencies * np.log(frequencies))
                    summary['shannon_diversity'] = shannon
                
                # Simpson diversity
                if total_reads > 0:
                    simpson = 1 - sum((df['cloneCount'] / total_reads) ** 2)
                    summary['simpson_diversity'] = simpson
                
                # Clonal expansion metrics
                if total_reads > 0:
                    # Top clone proportion
                    top_clone_count = df['cloneCount'].max()
                    summary['top_clone_proportion'] = top_clone_count / total_reads
                    
                    # Top 10 clones proportion
                    if total_clones >= 10:
                        top10_count = df.nlargest(10, 'cloneCount')['cloneCount'].sum()
                        summary['top10_proportion'] = top10_count / total_reads
                    else:
                        summary['top10_proportion'] = 1.0
            
            all_data.append(summary)
            print(f"Processed {patient}_{timepoint}: {total_clones} clones, {total_reads} reads")
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    if not all_data:
        print("No data could be processed")
        return None
    
    # Convert to DataFrame
    summary_df = pd.DataFrame(all_data)
    
    # Sort by patient and timepoint
    timepoint_order = {
        'T0_baseline': 0,
        'T1_cycle1': 1,
        'T2_cycle3': 2,
        'T3_progression': 3,
        'T4_posttreatment': 4
    }
    
    # Add timepoint order for sorting
    summary_df['timepoint_order'] = summary_df['timepoint'].map(
        lambda x: timepoint_order.get(x, 999)
    )
    
    # Sort by patient and timepoint
    summary_df = summary_df.sort_values(['patient', 'timepoint_order'])
    
    return summary_df

def create_diversity_plots(data, output_dir):
    """Create diversity plots."""
    if data is None or len(data) == 0:
        print("No data available for plotting")
        return
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with multiple subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('TCR Repertoire Analysis - 5 Timepoint Longitudinal Study', fontsize=16, fontweight='bold')
    
    # Plot 1: Shannon Diversity by Patient and Timepoint
    ax1 = axes[0, 0]
    for patient in data['patient'].unique():
        patient_data = data[data['patient'] == patient]
        if 'shannon_diversity' in patient_data.columns:
            ax1.plot(range(len(patient_data)), patient_data['shannon_diversity'], 'o-', 
                    label=patient, linewidth=2, markersize=8)
    ax1.set_title('Shannon Diversity Over Time', fontweight='bold')
    ax1.set_xlabel('Timepoint')
    ax1.set_ylabel('Shannon Diversity Index')
    ax1.set_xticks(range(len(data['timepoint'].unique())))
    ax1.set_xticklabels(data.sort_values('timepoint_order')['timepoint'].unique(), rotation=45)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Simpson Diversity
    ax2 = axes[0, 1]
    for patient in data['patient'].unique():
        patient_data = data[data['patient'] == patient]
        if 'simpson_diversity' in patient_data.columns:
            ax2.plot(range(len(patient_data)), patient_data['simpson_diversity'], 's-', 
                    label=patient, linewidth=2, markersize=8)
    ax2.set_title('Simpson Diversity Over Time', fontweight='bold')
    ax2.set_xlabel('Timepoint')
    ax2.set_ylabel('Simpson Diversity Index')
    ax2.set_xticks(range(len(data['timepoint'].unique())))
    ax2.set_xticklabels(data.sort_values('timepoint_order')['timepoint'].unique(), rotation=45)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Total Clones
    ax3 = axes[0, 2]
    for patient in data['patient'].unique():
        patient_data = data[data['patient'] == patient]
        ax3.plot(range(len(patient_data)), patient_data['total_clones'], '^-', 
                label=patient, linewidth=2, markersize=8)
    ax3.set_title('Total Clonotypes Over Time', fontweight='bold')
    ax3.set_xlabel('Timepoint')
    ax3.set_ylabel('Number of Clonotypes')
    ax3.set_xticks(range(len(data['timepoint'].unique())))
    ax3.set_xticklabels(data.sort_values('timepoint_order')['timepoint'].unique(), rotation=45)
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Top Clone Proportion
    ax4 = axes[1, 0]
    for patient in data['patient'].unique():
        patient_data = data[data['patient'] == patient]
        if 'top_clone_proportion' in patient_data.columns:
            ax4.plot(range(len(patient_data)), patient_data['top_clone_proportion'], 'o-', 
                    label=patient, linewidth=2, markersize=8)
    ax4.set_title('Top Clone Proportion Over Time', fontweight='bold')
    ax4.set_xlabel('Timepoint')
    ax4.set_ylabel('Top Clone Proportion')
    ax4.set_xticks(range(len(data['timepoint'].unique())))
    ax4.set_xticklabels(data.sort_values('timepoint_order')['timepoint'].unique(), rotation=45)
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: Top 10 Clone Proportion
    ax5 = axes[1, 1]
    for patient in data['patient'].unique():
        patient_data = data[data['patient'] == patient]
        if 'top10_proportion' in patient_data.columns:
            ax5.plot(range(len(patient_data)), patient_data['top10_proportion'], 's-', 
                    label=patient, linewidth=2, markersize=8)
    ax5.set_title('Top 10 Clones Proportion Over Time', fontweight='bold')
    ax5.set_xlabel('Timepoint')
    ax5.set_ylabel('Top 10 Clones Proportion')
    ax5.set_xticks(range(len(data['timepoint'].unique())))
    ax5.set_xticklabels(data.sort_values('timepoint_order')['timepoint'].unique(), rotation=45)
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: Diversity vs Clonal Expansion
    ax6 = axes[1, 2]
    if 'shannon_diversity' in data.columns and 'top_clone_proportion' in data.columns:
        scatter = ax6.scatter(data['shannon_diversity'], data['top_clone_proportion'], 
                            c=data['timepoint_order'], cmap='viridis', s=100, alpha=0.7)
        ax6.set_title('Diversity vs Clonal Expansion', fontweight='bold')
        ax6.set_xlabel('Shannon Diversity')
        ax6.set_ylabel('Top Clone Proportion')
        plt.colorbar(scatter, ax=ax6, label='Timepoint')
        ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'tcr_longitudinal_summary.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'tcr_longitudinal_summary.pdf'), bbox_inches='tight')
    
    print(f"✅ Plots saved to {output_dir}")
    
    # Save the data
    data.to_csv(os.path.join(output_dir, 'tcr_metrics.csv'), index=False)
    print(f"✅ Metrics saved to {os.path.join(output_dir, 'tcr_metrics.csv')}")

def create_simulated_data(output_dir):
    """Create simulated data for demonstration purposes."""
    print("⚠️ Creating simulated data for demonstration")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Define timepoints
    timepoints = ['T0_baseline', 'T1_cycle1', 'T2_cycle3', 'T3_progression', 'T4_posttreatment']
    timepoint_order = {tp: i for i, tp in enumerate(timepoints)}
    
    # Create simulated data
    data = []
    
    # Patient 1 - Expanding clones pattern
    for tp in timepoints:
        tp_idx = timepoint_order[tp]
        shannon = 4.5 - tp_idx * 0.5 + random.uniform(-0.2, 0.2)  # Decreasing diversity
        simpson = 0.95 - tp_idx * 0.05 + random.uniform(-0.02, 0.02)  # Decreasing diversity
        top_clone = 0.05 + tp_idx * 0.05 + random.uniform(-0.01, 0.01)  # Increasing dominance
        top10 = 0.2 + tp_idx * 0.1 + random.uniform(-0.02, 0.02)  # Increasing dominance
        
        data.append({
            'patient': 'PATIENT_01',
            'timepoint': tp,
            'timepoint_order': tp_idx,
            'sample': f'PATIENT_01_{tp}',
            'total_clones': int(1000 - tp_idx * 100 + random.uniform(-50, 50)),
            'total_reads': int(10000 + random.uniform(-500, 500)),
            'shannon_diversity': shannon,
            'simpson_diversity': simpson,
            'top_clone_proportion': top_clone,
            'top10_proportion': top10
        })
    
    # Patient 2 - Contracting clones pattern
    for tp in timepoints:
        tp_idx = timepoint_order[tp]
        shannon = 3.0 + tp_idx * 0.3 + random.uniform(-0.2, 0.2)  # Increasing diversity
        simpson = 0.7 + tp_idx * 0.05 + random.uniform(-0.02, 0.02)  # Increasing diversity
        top_clone = 0.25 - tp_idx * 0.03 + random.uniform(-0.01, 0.01)  # Decreasing dominance
        top10 = 0.6 - tp_idx * 0.08 + random.uniform(-0.02, 0.02)  # Decreasing dominance
        
        data.append({
            'patient': 'PATIENT_02',
            'timepoint': tp,
            'timepoint_order': tp_idx,
            'sample': f'PATIENT_02_{tp}',
            'total_clones': int(500 + tp_idx * 150 + random.uniform(-50, 50)),
            'total_reads': int(10000 + random.uniform(-500, 500)),
            'shannon_diversity': shannon,
            'simpson_diversity': simpson,
            'top_clone_proportion': top_clone,
            'top10_proportion': top10
        })
    
    # Convert to DataFrame
    df = pd.DataFrame(data)
    
    # Sort by patient and timepoint
    df = df.sort_values(['patient', 'timepoint_order'])
    
    return df

def main():
    """Main function."""
    print("🧬 TCR ANALYSIS PLOT GENERATOR")
    print("==============================")
    
    # Define directories
    input_dir = "tcr_analysis_results"
    output_dir = "tcr_analysis_results/plots"
    
    # Check if input directory exists
    if not os.path.exists(input_dir):
        print(f"⚠️ Input directory {input_dir} not found")
        print("Creating simulated data for demonstration")
        data = create_simulated_data(output_dir)
    else:
        # Load data
        data = load_clonotype_data(input_dir)
        
        # If no data, create simulated data
        if data is None or len(data) == 0:
            print("⚠️ No data found in input directory")
            print("Creating simulated data for demonstration")
            data = create_simulated_data(output_dir)
    
    # Create plots
    create_diversity_plots(data, output_dir)
    
    print("\n🎉 PLOT GENERATION COMPLETE!")
    print("============================")
    print(f"📁 Plots saved to: {output_dir}")
    print("📊 Files generated:")
    print(f"   • {output_dir}/tcr_longitudinal_summary.png")
    print(f"   • {output_dir}/tcr_longitudinal_summary.pdf")
    print(f"   • {output_dir}/tcr_metrics.csv")
    print("\n🔍 To view plots:")
    print(f"   open {output_dir}/tcr_longitudinal_summary.pdf")

if __name__ == "__main__":
    main()
