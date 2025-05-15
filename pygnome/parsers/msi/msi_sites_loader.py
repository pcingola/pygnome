"""
Utility functions for loading MSI sites into a GenomicFeatureStore.
"""

from pygnome.feature_store.genomic_feature_store import GenomicFeatureStore, StoreType
from pygnome.feature_store.msi_chromosome_store import MsiChromosomeStore, MsiSiteCounter
from pygnome.parsers.msi.msi_sites_reader import MsiSitesReader


def load_msi_sites(file_path: str) -> GenomicFeatureStore:
    """
    Load MSI sites from a file into a GenomicFeatureStore.
    
    Uses a two-pass approach for memory efficiency:
    1. First pass: Count features per chromosome and calculate statistics
    2. Second pass: Load features into pre-allocated arrays
    
    Args:
        file_path: Path to the MSI sites file
        
    Returns:
        A GenomicFeatureStore containing the MSI sites
    """
    # Create the genomic feature store
    feature_store = GenomicFeatureStore()
    
    # First pass: Count features per chromosome
    counter = MsiSiteCounter()
    with MsiSitesReader(file_path) as reader:
        for record in reader:
            counter.add(record)
    
    # Create chromosome stores with pre-allocated arrays
    for chrom in counter.feature_counts:
        count = counter.get_count(chrom)
        max_lengths = counter.get_max_lengths(chrom)
        
        # Create and register the chromosome store
        chrom_store = MsiChromosomeStore(feature_count=count, max_lengths_by_bin=max_lengths)
        feature_store.chromosomes[chrom] = chrom_store
        chrom_store.index_build_start()
    
    # Second pass: Load features into arrays
    with MsiSitesReader(file_path) as reader:
        for record in reader:
            chrom_store = feature_store.chromosomes.get(record.chrom)
            if chrom_store:
                chrom_store.add(record)
    
    # Finalize all chromosome stores
    for chrom_store in feature_store.chromosomes.values():
        chrom_store.index_build_end()
    
    return feature_store