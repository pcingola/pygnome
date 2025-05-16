#!/usr/bin/env python3
"""
Diagnostic script to analyze why the pickled GenomicFeatureStore is larger than expected.
Run this script after creating your GenomicFeatureStore but before saving it.
"""

import sys
import os
import numpy as np
import pickle
from pathlib import Path

from pygnome.feature_store.msi_chromosome_store import MsiChromosomeStore
from pygnome.parsers.msi.msi_sites_loader import load_msi_sites

def analyze_store(gfs):
    """Analyze the memory usage of a GenomicFeatureStore."""
    print("\n--- Memory Usage Analysis ---")
    total_allocated = 0
    total_used = 0
    
    for chrom, store in gfs.chromosomes.items():
        if isinstance(store, MsiChromosomeStore):
            # Calculate allocated space
            starts_size = store._starts.nbytes
            ends_size = store._ends.nbytes
            
            # Calculate used space
            used_starts = store._feature_count * np.dtype(np.uint32).itemsize
            used_ends = store._feature_count * np.dtype(np.uint32).itemsize
            
            # Get DnaStringArray stats
            dna_count, dna_total_nt, dna_bits_per_nt = store._repeat_unit_bases.get_stats()
            dna_allocated = (len(store._repeat_unit_bases._data) * 4) + (len(store._repeat_unit_bases._positions) * 8)
            
            print(f"Chromosome {chrom}:")
            print(f"  Features: {store._feature_count}")
            efficiency_starts = used_starts/starts_size if starts_size > 0 else 0
            efficiency_ends = used_ends/ends_size if ends_size > 0 else 0
            print(f"  Starts array: {starts_size} bytes allocated, {used_starts} bytes used ({efficiency_starts:.2%})")
            print(f"  Ends array: {ends_size} bytes allocated, {used_ends} bytes used ({efficiency_ends:.2%})")
            print(f"  DnaStringArray: {dna_allocated} bytes allocated")
            
            total_allocated += starts_size + ends_size + dna_allocated
            total_used += used_starts + used_ends

    print(f"Total allocated: {total_allocated} bytes")
    print(f"Total used: {total_used} bytes")
    efficiency = total_used/total_allocated if total_allocated > 0 else 0
    print(f"Efficiency: {efficiency:.2%}")
    
    return total_allocated, total_used

def trim_arrays(gfs):
    """Trim arrays in MsiChromosomeStore to their actual size."""
    print("\n--- Trimming Arrays ---")
    for chrom, store in gfs.chromosomes.items():
        if isinstance(store, MsiChromosomeStore):
            if store._feature_count < len(store._starts):
                print(f"Trimming arrays for {chrom}...")
                # Trim the arrays
                store._starts = store._starts[:store._feature_count].copy()
                store._ends = store._ends[:store._feature_count].copy()
                
                # We can't easily trim the DnaStringArray, but we could create a new one
                # This would require more complex code to copy all the strings
    
    return gfs

def compare_file_sizes(original_file, pickle_file):
    """Compare the sizes of the original and pickled files."""
    original_size = os.path.getsize(original_file)
    pickle_size = os.path.getsize(pickle_file)
    
    print(f"\n--- File Size Comparison ---")
    print(f"Original file size: {original_size} bytes")
    print(f"Pickled file size: {pickle_size} bytes")
    print(f"Ratio: {pickle_size/original_size:.2f}x larger")
    
    return original_size, pickle_size

def main():
    """Main function to run the analysis."""
    if len(sys.argv) < 2:
        print("Usage: python analyze_msi_store.py <msi_sites_file>")
        sys.exit(1)
    
    msi_sites_file = Path(sys.argv[1])
    save_to = msi_sites_file.with_suffix('.pckl')
    trimmed_save_to = msi_sites_file.with_suffix('.trimmed.pckl')
    
    # Load the MSI sites
    print(f"Loading MSI sites from {msi_sites_file}...")
    gfs = load_msi_sites(msi_sites_file, verbose=True)
    
    # Analyze the store
    print("\nAnalyzing original store...")
    analyze_store(gfs)
    
    # Save the original store
    print(f"\nSaving original store to {save_to}...")
    gfs.save(save_to)
    
    # Trim the arrays
    print("\nTrimming arrays...")
    trimmed_gfs = trim_arrays(gfs)
    
    # Analyze the trimmed store
    print("\nAnalyzing trimmed store...")
    analyze_store(trimmed_gfs)
    
    # Save the trimmed store
    print(f"\nSaving trimmed store to {trimmed_save_to}...")
    trimmed_gfs.save(trimmed_save_to)
    
    # Compare file sizes
    print("\nComparing original file sizes...")
    compare_file_sizes(msi_sites_file, save_to)
    
    print("\nComparing trimmed file sizes...")
    compare_file_sizes(msi_sites_file, trimmed_save_to)

if __name__ == "__main__":
    main()