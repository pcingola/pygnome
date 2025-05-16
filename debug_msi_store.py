#!/usr/bin/env python3
"""
Debug script for investigating the MSI store issue.
Uses a tiny MSI sites file for easier debugging.
"""

from pathlib import Path
import os
import sys
import numpy as np

from pygnome.feature_store.genomic_feature_store import GenomicFeatureStore, StoreType
from pygnome.feature_store.msi_chromosome_store import MsiChromosomeStore, MsiSiteCounter
from pygnome.parsers.msi.msi_sites_reader import MsiSitesReader
from pygnome.parsers.msi.msi_sites_loader import load_msi_sites

def debug_load_msi_sites(file_path: Path, verbose: bool = True) -> GenomicFeatureStore:
    """
    Debug version of load_msi_sites with extra logging.
    """
    # Create the genomic feature store
    feature_store = GenomicFeatureStore(store_type=StoreType.MSI)
    
    # First pass: Count features per chromosome
    counter = MsiSiteCounter()
    if verbose:
        print(f"Counting features in {file_path}...")
    
    # Debug: Print the first few records
    print("\n--- First Pass: Counting Features ---")
    count = 0
    with MsiSitesReader(file_path) as reader:
        for record in reader:
            print(f"Record {count}: chrom={record.chrom}, start={record.start}, end={record.end}")
            print(f"  type={type(record).__name__}, repeat_unit_bases={record.repeat_unit_bases}")
            counter.add(record)
            count += 1
    
    print(f"Total features counted: {count}")
    
    # Create chromosome stores with pre-allocated arrays
    print("\n--- Creating Chromosome Stores ---")
    for chrom in counter.feature_counts:
        count = counter.get_count(chrom)
        max_lengths = counter.get_max_lengths(chrom)
        
        print(f"Creating store for {chrom}: {count} features")
        
        # Create and register the chromosome store
        chrom_store = MsiChromosomeStore(chrom=chrom, feature_count=count, max_lengths_by_bin=max_lengths)
        feature_store.chromosomes[chrom] = chrom_store
        chrom_store.index_build_start()
    
    # Second pass: Load features into arrays
    print("\n--- Second Pass: Loading Features ---")
    count = 0
    # Print the keys in the chromosomes dictionary
    print("\n--- Chromosome Store Keys ---")
    print(f"Dictionary: {feature_store.chromosomes}")
    for key in feature_store.chromosomes.keys():
        print(f"  Key: '{key}' (type: {type(key).__name__}, id: {id(key)})")
    
    with MsiSitesReader(file_path) as reader:
        for record in reader:
            chrom_key = record.chrom
            print(f"Adding record {count}: chrom='{chrom_key}'")
            # Add the record to the store
            feature_store.add(record)
            count += 1
    
    print(f"Total features loaded: {count}")
    
    # Check final counts
    print("\n--- Final Feature Counts ---")
    total_features = 0
    for chrom, store in feature_store.chromosomes.items():
        print(f"  {chrom}: {len(store)} features")
        total_features += len(store)
    print(f"  Total features in stores: {total_features}")
    
    # Finalize all chromosome stores
    for chrom_store in feature_store.chromosomes.values():
        chrom_store.index_build_end()
    
    return feature_store

def main():
    """Main function to run the debug script."""
    # Use the tiny file for debugging
    home_path = Path.home()
    msi_sites_file = home_path / 'hg38_sites.tiny.txt'
    
    if not msi_sites_file.exists():
        print(f"Error: File not found: {msi_sites_file}")
        sys.exit(1)
    
    print(f"Using MSI sites file: {msi_sites_file}")
    
    # Load the MSI sites with debug logging
    gfs = debug_load_msi_sites(msi_sites_file)
    
    # Analyze memory usage before saving
    print("\n--- Memory Analysis Before Saving ---")
    for chrom, store in gfs.chromosomes.items():
        if isinstance(store, MsiChromosomeStore):
            starts_allocated = store._starts.nbytes
            ends_allocated = store._ends.nbytes
            starts_used = store._feature_count * np.dtype(np.uint32).itemsize
            ends_used = store._feature_count * np.dtype(np.uint32).itemsize
            
            print(f"Chromosome {chrom}:")
            print(f"  Features: {store._feature_count}")
            print(f"  Starts array: {starts_allocated} bytes allocated, {starts_used} bytes used ({starts_used/starts_allocated:.2%})")
            print(f"  Ends array: {ends_allocated} bytes allocated, {ends_used} bytes used ({ends_used/ends_allocated:.2%})")
            
            if hasattr(store, '_repeat_unit_bases'):
                dna_array = store._repeat_unit_bases
                dna_count, dna_total_nt, dna_bits_per_nt = dna_array.get_stats()
                data_bytes = len(dna_array._data) * 4  # 4 bytes per uint32
                positions_bytes = len(dna_array._positions) * 8  # 8 bytes per uint64
                print(f"  DnaStringArray:")
                print(f"    Data array: {data_bytes} bytes allocated")
                print(f"    Positions array: {positions_bytes} bytes allocated")
                print(f"    Total nucleotides: {dna_total_nt}")
                print(f"    Bits per nucleotide: {dna_bits_per_nt:.2f}")
    
    # Save the store
    save_to = msi_sites_file.with_suffix('.pckl')
    print(f"\nSaving store to {save_to}...")
    gfs.save(save_to)
    
    # Compare file sizes
    original_size = os.path.getsize(msi_sites_file)
    saved_size = os.path.getsize(save_to)
    
    print("\n--- File Size Comparison ---")
    print(f"Original file size: {original_size} bytes")
    print(f"Saved file size: {saved_size} bytes")
    print(f"Ratio: {saved_size/original_size:.2f}x larger")
    
    # Create a manually trimmed version for comparison
    print("\n--- Creating Manually Trimmed Version ---")
    # First, let's analyze the DnaStringArray in detail
    for chrom, store in gfs.chromosomes.items():
        if isinstance(store, MsiChromosomeStore) and hasattr(store, '_repeat_unit_bases'):
            dna_array = store._repeat_unit_bases
            data_bytes = len(dna_array._data) * 4  # 4 bytes per uint32
            positions_bytes = len(dna_array._positions) * 8  # 8 bytes per uint64
            print(f"DnaStringArray for {chrom} before trim:")
            print(f"  Data array: {len(dna_array._data)} elements, {data_bytes} bytes")
            print(f"  Positions array: {len(dna_array._positions)} elements, {positions_bytes} bytes")
            print(f"  Total nucleotides: {dna_array._total_nucleotides}")
            
            # Now call the trim method
            print("Calling trim() method...")
            store.trim()
            
            # Check the size after trimming
            data_bytes_after = len(dna_array._data) * 4
            positions_bytes_after = len(dna_array._positions) * 8
            print(f"DnaStringArray for {chrom} after trim:")
            print(f"  Data array: {len(dna_array._data)} elements, {data_bytes_after} bytes")
            print(f"  Positions array: {len(dna_array._positions)} elements, {positions_bytes_after} bytes")
            print(f"  Memory reduction: {(data_bytes + positions_bytes - data_bytes_after - positions_bytes_after) / (data_bytes + positions_bytes):.2%}")
    
    # Save the trimmed store
    trimmed_save_to = msi_sites_file.with_suffix('.trimmed.pckl')
    print(f"\nSaving trimmed store to {trimmed_save_to}...")
    gfs.save(trimmed_save_to)
    
    # Compare with trimmed version
    if os.path.exists(trimmed_save_to):
        trimmed_size = os.path.getsize(trimmed_save_to)
        print(f"Trimmed file size: {trimmed_size} bytes")
        print(f"Trimmed ratio: {trimmed_size/original_size:.2f}x original")
        print(f"Improvement: {(saved_size - trimmed_size) / saved_size:.2%} smaller than untrimmed")
        
    # Create a new store with smaller initial allocation
    print("\n--- Creating Store with Smaller Initial Allocation ---")
    # Monkey patch the DEFAULT_CAPACITY constant to a smaller value
    import pygnome.sequences.dna_string_array
    original_capacity = pygnome.sequences.dna_string_array.DEFAULT_CAPACITY
    pygnome.sequences.dna_string_array.DEFAULT_CAPACITY = 1024  # Use 1KB instead of 1MB
    
    # Load the data again
    print("Loading data with smaller allocation...")
    small_gfs = debug_load_msi_sites(msi_sites_file)
    
    # Save the small allocation store
    small_save_to = msi_sites_file.with_suffix('.small.pckl')
    print(f"Saving small allocation store to {small_save_to}...")
    small_gfs.save(small_save_to)
    
    # Compare all versions
    if os.path.exists(small_save_to):
        small_size = os.path.getsize(small_save_to)
        print("\n--- Final Comparison ---")
        print(f"Original text file: {original_size} bytes")
        print(f"Default pickle file: {saved_size} bytes ({saved_size/original_size:.2f}x)")
        print(f"Trimmed pickle file: {trimmed_size} bytes ({trimmed_size/original_size:.2f}x)")
        print(f"Small allocation pickle file: {small_size} bytes ({small_size/original_size:.2f}x)")
        
    # Restore the original capacity
    pygnome.sequences.dna_string_array.DEFAULT_CAPACITY = original_capacity

if __name__ == "__main__":
    main()