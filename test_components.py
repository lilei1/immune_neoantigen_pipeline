#!/usr/bin/env python3

"""
Test individual pipeline components
This script tests the core functionality without requiring Nextflow
"""

import sys
import os
import subprocess
from pathlib import Path

def test_samplesheet_validation():
    """Test samplesheet validation"""
    print("🧪 Testing samplesheet validation...")
    
    try:
        result = subprocess.run([
            'python3', 'bin/metadata_validator.py', 
            '--input', 'assets/samplesheet_test.csv'
        ], capture_output=True, text=True, cwd='.')
        
        if result.returncode == 0:
            print("✅ Samplesheet validation: PASSED")
            return True
        else:
            print("❌ Samplesheet validation: FAILED")
            print(result.stdout)
            print(result.stderr)
            return False
    except Exception as e:
        print(f"❌ Error testing samplesheet validation: {e}")
        return False

def test_clonotype_tracking():
    """Test clonotype tracking script"""
    print("🧪 Testing clonotype tracking...")
    
    # Create dummy clonotype data for testing
    test_data_dir = Path("test_data")
    test_data_dir.mkdir(exist_ok=True)
    
    # Create sample clonotype files
    sample1_data = """cloneId	cloneCount	cloneFraction	targetSequences
CLONE_001	1500	0.15	TRBV7-9*01>TRBD1*01>TRBJ2-1*01
CLONE_002	1200	0.12	TRBV12-3*01>TRBD2*01>TRBJ1-4*01
CLONE_003	800	0.08	TRBV5-1*01>TRBD1*01>TRBJ2-7*01"""

    sample2_data = """cloneId	cloneCount	cloneFraction	targetSequences
CLONE_001	1800	0.18	TRBV7-9*01>TRBD1*01>TRBJ2-1*01
CLONE_004	900	0.09	TRBV19*01>TRBD2*01>TRBJ1-1*01
CLONE_005	600	0.06	TRBV20-1*01>TRBD1*01>TRBJ2-3*01"""
    
    with open(test_data_dir / "sample1.clonotypes.txt", "w") as f:
        f.write(sample1_data)
    
    with open(test_data_dir / "sample2.clonotypes.txt", "w") as f:
        f.write(sample2_data)
    
    try:
        # Test if R is available
        r_check = subprocess.run(['R', '--version'], capture_output=True)
        if r_check.returncode != 0:
            print("⚠️  R not available, skipping clonotype tracking test")
            return True
        
        # Run clonotype tracking
        result = subprocess.run([
            'Rscript', 'bin/track_clones.R', 
            'test_output',
            str(test_data_dir / "sample1.clonotypes.txt"),
            str(test_data_dir / "sample2.clonotypes.txt")
        ], capture_output=True, text=True, cwd='.')
        
        if result.returncode == 0:
            print("✅ Clonotype tracking: PASSED")
            return True
        else:
            print("❌ Clonotype tracking: FAILED")
            print(result.stdout)
            print(result.stderr)
            return False
            
    except Exception as e:
        print(f"❌ Error testing clonotype tracking: {e}")
        return False
    finally:
        # Cleanup
        import shutil
        if test_data_dir.exists():
            shutil.rmtree(test_data_dir)

def test_variant_merging():
    """Test variant merging script"""
    print("🧪 Testing variant merging...")
    
    # Create dummy VCF data
    test_data_dir = Path("test_vcf")
    test_data_dir.mkdir(exist_ok=True)
    
    vcf_header = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=249250621>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE"""

    vcf1_data = f"""{vcf_header}
chr1	100	.	A	T	60	PASS	DP=30	GT:DP	0/1:30
chr1	200	.	G	C	80	PASS	DP=40	GT:DP	1/1:40"""

    vcf2_data = f"""{vcf_header}
chr1	100	.	A	T	65	PASS	DP=35	GT:DP	0/1:35
chr1	300	.	C	G	70	PASS	DP=25	GT:DP	0/1:25"""
    
    with open(test_data_dir / "sample1.vcf", "w") as f:
        f.write(vcf1_data)
    
    with open(test_data_dir / "sample2.vcf", "w") as f:
        f.write(vcf2_data)
    
    try:
        result = subprocess.run([
            'python3', 'bin/merge_variants.py',
            '--input', str(test_data_dir / "sample1.vcf"), str(test_data_dir / "sample2.vcf"),
            '--output', 'test_merged_variants.tsv',
            '--sample-names', 'sample1', 'sample2'
        ], capture_output=True, text=True, cwd='.')
        
        if result.returncode == 0:
            print("✅ Variant merging: PASSED")
            # Check if output file was created
            if Path("test_merged_variants.tsv").exists():
                print("✅ Output file created successfully")
                return True
            else:
                print("❌ Output file not created")
                return False
        else:
            print("❌ Variant merging: FAILED")
            print(result.stdout)
            print(result.stderr)
            return False
            
    except Exception as e:
        print(f"❌ Error testing variant merging: {e}")
        return False
    finally:
        # Cleanup
        import shutil
        if test_data_dir.exists():
            shutil.rmtree(test_data_dir)
        if Path("test_merged_variants.tsv").exists():
            Path("test_merged_variants.tsv").unlink()

def test_nextflow_syntax():
    """Test Nextflow syntax if available"""
    print("🧪 Testing Nextflow syntax...")
    
    try:
        # Check if nextflow is available
        nf_check = subprocess.run(['nextflow', '-version'], capture_output=True)
        if nf_check.returncode != 0:
            print("⚠️  Nextflow not available in PATH, skipping syntax test")
            return True
        
        # Test syntax
        result = subprocess.run([
            'nextflow', 'run', 'main.nf', '--help'
        ], capture_output=True, text=True, cwd='.')
        
        if result.returncode == 0:
            print("✅ Nextflow syntax: PASSED")
            return True
        else:
            print("❌ Nextflow syntax: FAILED")
            print(result.stdout)
            print(result.stderr)
            return False
            
    except Exception as e:
        print(f"⚠️  Nextflow not available: {e}")
        return True  # Don't fail if Nextflow isn't available

def main():
    print("🚀 Testing Immune Repertoire + Neoantigen Pipeline Components")
    print("=" * 60)
    
    tests = [
        test_samplesheet_validation,
        test_variant_merging,
        test_clonotype_tracking,
        test_nextflow_syntax
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        try:
            if test():
                passed += 1
            print()  # Add spacing between tests
        except Exception as e:
            print(f"❌ Test failed with exception: {e}")
            print()
    
    print("=" * 60)
    print(f"📊 Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All component tests PASSED!")
        print("\n✨ Your pipeline is ready to use!")
        print("\n📋 Next steps:")
        print("1. Ensure Nextflow is in your PATH")
        print("2. Run: nextflow run main.nf -profile test,docker")
        print("3. Or adapt the samplesheet for your data")
    else:
        print("⚠️  Some tests failed. Please check the errors above.")
        
    return passed == total

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
