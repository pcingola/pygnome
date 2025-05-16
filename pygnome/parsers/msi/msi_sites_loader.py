"""
Utility functions for loading MSI sites into a GenomicFeatureStore.
"""

from pathlib import Path
import time
from pygnome.feature_store.genomic_feature_store import GenomicFeatureStore, StoreType
from pygnome.feature_store.msi_chromosome_store import MsiChromosomeStore, MsiSiteCounter
from pygnome.parsers.msi.msi_sites_reader import MsiSitesReader


def load_msi_sites(file_path: Path, verbose: bool, force: bool = False) -> GenomicFeatureStore:
    """
    Load MSI sites from a file into a GenomicFeatureStore.
    
    Uses a two-pass approach for memory efficiency:
    1. First pass: Count features per chromosome and calculate statistics
    2. Second pass: Load features into pre-allocated arrays
    
    Args:
        file_path: Path to the MSI sites file
        verbose: If True, print progress messages

    Returns:
        A GenomicFeatureStore containing the MSI sites
    """
    # Create the genomic feature store
    feature_store = GenomicFeatureStore(store_type=StoreType.MSI)

    pickle_path = file_path.with_suffix('.pckl')
    if pickle_path.exists() and not force:
        # Load from the pickle file if it exists and force is not set
        if verbose:
            print(f"Loading from {pickle_path}...")
        return GenomicFeatureStore.load(pickle_path)

    # First pass: Count features per chromosome
    counter = MsiSiteCounter()
    if verbose:
        print(f"Counting features in {file_path}...")
    prev_chrom, count = None, 0
    start_time = time.time()
    with MsiSitesReader(file_path) as reader:
        for record in reader:
            counter.add(record)
            count += 1
            if verbose and record.chrom != prev_chrom:
                # Show only if new chromosome is encountered
                elapsed = time.time() - start_time
                print(f"{count:,} ({elapsed:.2f}s): {record.chrom}")
                prev_chrom = record.chrom
    if verbose:
        elapsed = time.time() - start_time
        print(f"Total features counted: {count} (elapsed: {elapsed:.2f}s)")

    # Create chromosome stores with pre-allocated arrays
    for chrom in counter.feature_counts:
        count = counter.get_count(chrom)
        max_lengths = counter.get_max_lengths(chrom)
        
        # Create and register the chromosome store
        chrom_store = MsiChromosomeStore(chrom=chrom, feature_count=count, max_lengths_by_bin=max_lengths)
        feature_store.chromosomes[chrom] = chrom_store
        chrom_store.index_build_start()
    
    # Second pass: Load features into arrays
    if verbose:
        print(f"Loading features into stores...")
    prev_chrom, count = None, 0
    with MsiSitesReader(file_path) as reader:
        for record in reader:
            # Use feature_store.add() instead of directly accessing the chromosome store
            feature_store.add(record)
            count += 1
            if verbose and record.chrom != prev_chrom:
                # Show only if new chromosome is encountered
                elapsed = time.time() - start_time
                print(f"{count:,} ({elapsed:.2f}s): {record.chrom}")
                prev_chrom = record.chrom
    if verbose:
        elapsed = time.time() - start_time
        print(f"Total features loaded: {count} ({elapsed:.2f}s)")
        
        # Debug: Check if features were actually added
        total_features = 0
        for chrom, store in feature_store.chromosomes.items():
            print(f"  {chrom}: {len(store)} features")
            total_features += len(store)
        print(f"  Total features in store: {total_features}")

    # Finalize all chromosome stores
    for chrom_store in feature_store.chromosomes.values():
        if verbose:
            print(f"Indexing '{chrom_store.chromosome}'...")
        chrom_store.index_build_end()
    
    # Save the feature store to a pickle file
    if verbose:
        print(f"Saving to {pickle_path}...")
    feature_store.save(pickle_path)
    if verbose:
        print(f"Saved to {pickle_path}, file size: {pickle_path.stat().st_size / 1e6:.2f} MB")

    return feature_store