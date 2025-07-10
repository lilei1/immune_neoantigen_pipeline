#!/usr/bin/env python3
"""
Demonstrate TCR filtering and collapsing workflow with MiXCR.
This script shows the complete pipeline from raw reads to filtered/collapsed clonotypes.
"""

import subprocess
import os
import pandas as pd
from pathlib import Path

def run_command(cmd, description):
    """Run a command and handle errors."""
    print(f"\n🔄 {description}")
    print(f"Command: {cmd}")
    
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        print(f"✅ Success: {description}")
        if result.stdout:
            print(f"Output: {result.stdout[:200]}...")
        return result
    except subprocess.CalledProcessError as e:
        print(f"❌ Error: {description}")
        print(f"Error output: {e.stderr}")
        return None

def analyze_clonotype_file(file_path, description):
    """Analyze a clonotype file and print statistics."""
    if not os.path.exists(file_path):
        print(f"❌ File not found: {file_path}")
        return None
    
    try:
        # Try to read the file
        if file_path.endswith('.tsv'):
            df = pd.read_csv(file_path, sep='\t')
        else:
            df = pd.read_csv(file_path, sep='\t')  # MiXCR output is usually tab-separated
        
        print(f"\n📊 {description} - {file_path}")
        print(f"   Total clonotypes: {len(df)}")
        
        if 'cloneCount' in df.columns:
            print(f"   Total reads: {df['cloneCount'].sum():,}")
            print(f"   Mean clone size: {df['cloneCount'].mean():.1f}")
            print(f"   Median clone size: {df['cloneCount'].median():.1f}")
            print(f"   Top 10 clones represent: {df['cloneCount'].nlargest(10).sum() / df['cloneCount'].sum() * 100:.1f}% of reads")
        
        if 'cloneFraction' in df.columns:
            print(f"   Top clone frequency: {df['cloneFraction'].max():.4f}")
            print(f"   Clones >1% frequency: {(df['cloneFraction'] > 0.01).sum()}")
            print(f"   Clones >0.1% frequency: {(df['cloneFraction'] > 0.001).sum()}")
        
        return df
    except Exception as e:
        print(f"❌ Error reading {file_path}: {e}")
        return None

def demonstrate_tcr_workflow():
    """Demonstrate the complete TCR filtering and collapsing workflow."""
    
    print("🧬 TCR FILTERING AND COLLAPSING DEMONSTRATION")
    print("=" * 60)
    
    # Check if test data exists
    test_data_dir = Path("test_data_tcr")
    if not test_data_dir.exists():
        print("❌ Test data not found. Please run: python3 scripts/generate_tcr_longitudinal_data.py")
        return
    
    # Create output directory
    output_dir = Path("tcr_demo_output")
    output_dir.mkdir(exist_ok=True)
    
    # Select one sample for demonstration
    sample_name = "PATIENT_01_T0_baseline"
    r1_file = test_data_dir / f"{sample_name}_tcr_1.fastq.gz"
    r2_file = test_data_dir / f"{sample_name}_tcr_2.fastq.gz"
    
    if not r1_file.exists() or not r2_file.exists():
        print(f"❌ Sample files not found: {r1_file}, {r2_file}")
        return
    
    print(f"📁 Working with sample: {sample_name}")
    print(f"📁 Output directory: {output_dir}")
    
    # Step 1: MiXCR Align
    vdjca_file = output_dir / f"{sample_name}.vdjca"
    align_cmd = f"""docker run --rm -v $(pwd):/data -w /data mgibio/mixcr:latest \\
        mixcr align -s hsa \\
        --report /data/{output_dir}/{sample_name}_align.report \\
        /data/{r1_file} /data/{r2_file} \\
        /data/{vdjca_file}"""
    
    result = run_command(align_cmd, "MiXCR Alignment")
    if not result:
        return
    
    # Step 2: MiXCR Assemble
    clns_file = output_dir / f"{sample_name}.clns"
    assemble_cmd = f"""docker run --rm -v $(pwd):/data -w /data mgibio/mixcr:latest \\
        mixcr assemble \\
        --report /data/{output_dir}/{sample_name}_assemble.report \\
        /data/{vdjca_file} \\
        /data/{clns_file}"""
    
    result = run_command(assemble_cmd, "MiXCR Assembly")
    if not result:
        return
    
    # Step 3: Export raw clonotypes (before filtering)
    raw_clonotypes = output_dir / f"{sample_name}_raw.clonotypes.tsv"
    export_raw_cmd = f"""docker run --rm -v $(pwd):/data -w /data mgibio/mixcr:latest \\
        mixcr exportClones --chains ALL --preset full -f \\
        /data/{clns_file} \\
        /data/{raw_clonotypes}"""
    
    result = run_command(export_raw_cmd, "Export Raw Clonotypes")
    if not result:
        return
    
    # Analyze raw clonotypes
    raw_df = analyze_clonotype_file(raw_clonotypes, "RAW CLONOTYPES (Before Filtering)")
    
    # Step 4: MiXCR Filter
    filtered_clns = output_dir / f"{sample_name}_filtered.clns"
    filter_cmd = f"""docker run --rm -v $(pwd):/data -w /data mgibio/mixcr:latest \\
        mixcr filter \\
        --min-count 2 \\
        --min-fraction 0.00001 \\
        --report /data/{output_dir}/{sample_name}_filter.report \\
        /data/{clns_file} \\
        /data/{filtered_clns}"""
    
    result = run_command(filter_cmd, "MiXCR Filtering")
    if not result:
        return
    
    # Step 5: Export filtered clonotypes
    filtered_clonotypes = output_dir / f"{sample_name}_filtered.clonotypes.tsv"
    export_filtered_cmd = f"""docker run --rm -v $(pwd):/data -w /data mgibio/mixcr:latest \\
        mixcr exportClones --chains ALL --preset full -f \\
        /data/{filtered_clns} \\
        /data/{filtered_clonotypes}"""
    
    result = run_command(export_filtered_cmd, "Export Filtered Clonotypes")
    if not result:
        return
    
    # Analyze filtered clonotypes
    filtered_df = analyze_clonotype_file(filtered_clonotypes, "FILTERED CLONOTYPES (After Filtering)")
    
    # Step 6: MiXCR Collapse (Note: This command might not work with all MiXCR versions)
    collapsed_clns = output_dir / f"{sample_name}_collapsed.clns"
    collapse_cmd = f"""docker run --rm -v $(pwd):/data -w /data mgibio/mixcr:latest \\
        mixcr collapse \\
        --by-CDR3 \\
        --level 0.9 \\
        --report /data/{output_dir}/{sample_name}_collapse.report \\
        /data/{filtered_clns} \\
        /data/{collapsed_clns}"""
    
    result = run_command(collapse_cmd, "MiXCR Collapsing")
    
    if result:
        # Step 7: Export collapsed clonotypes
        collapsed_clonotypes = output_dir / f"{sample_name}_collapsed.clonotypes.tsv"
        export_collapsed_cmd = f"""docker run --rm -v $(pwd):/data -w /data mgibio/mixcr:latest \\
            mixcr exportClones --chains ALL --preset full -f \\
            /data/{collapsed_clns} \\
            /data/{collapsed_clonotypes}"""
        
        result = run_command(export_collapsed_cmd, "Export Collapsed Clonotypes")
        if result:
            # Analyze collapsed clonotypes
            collapsed_df = analyze_clonotype_file(collapsed_clonotypes, "COLLAPSED CLONOTYPES (After Collapsing)")
    else:
        print("⚠️  Collapsing step failed - this may not be available in this MiXCR version")
        print("   Filtering alone provides significant quality improvement")
    
    # Summary comparison
    print("\n" + "=" * 60)
    print("📈 FILTERING AND COLLAPSING SUMMARY")
    print("=" * 60)
    
    if raw_df is not None and filtered_df is not None:
        raw_count = len(raw_df)
        filtered_count = len(filtered_df)
        reduction = (raw_count - filtered_count) / raw_count * 100
        
        print(f"Raw clonotypes: {raw_count:,}")
        print(f"Filtered clonotypes: {filtered_count:,}")
        print(f"Reduction by filtering: {reduction:.1f}%")
        
        if 'cloneCount' in raw_df.columns and 'cloneCount' in filtered_df.columns:
            raw_reads = raw_df['cloneCount'].sum()
            filtered_reads = filtered_df['cloneCount'].sum()
            read_retention = filtered_reads / raw_reads * 100
            print(f"Read retention after filtering: {read_retention:.1f}%")
    
    print(f"\n📁 All output files saved to: {output_dir}/")
    print("🔍 Check the *_report files for detailed MiXCR statistics")
    
    # List output files
    print(f"\n📋 Generated files:")
    for file in sorted(output_dir.glob("*")):
        size = file.stat().st_size
        if size > 1024*1024:
            size_str = f"{size/(1024*1024):.1f}MB"
        elif size > 1024:
            size_str = f"{size/1024:.1f}KB"
        else:
            size_str = f"{size}B"
        print(f"   {file.name} ({size_str})")

if __name__ == "__main__":
    demonstrate_tcr_workflow()
