"""Tests for serialization of GenomicFeatureStore."""

import unittest
import tempfile
from pathlib import Path
import time

from pygnome.feature_store.genomic_feature_store import GenomicFeatureStore, StoreType
from pygnome.genomics.genomic_feature import GenomicFeature
from pygnome.genomics.strand import Strand
from pygnome.parsers.msi.msi_site_record import MsiSiteRecord


class TestFeatureStoreSerialization(unittest.TestCase):
    """Test cases for serialization of GenomicFeatureStore."""

    def setUp(self):
        """Set up test fixtures."""
        # Create a sample feature store
        self.store = GenomicFeatureStore(store_type=StoreType.INTERVAL_TREE)
        
        # Add some features
        with self.store:
            # Chromosome 1
            self.store.add(GenomicFeature(id="gene1", chrom="chr1", start=100, end=200, strand=Strand.POSITIVE))
            self.store.add(GenomicFeature(id="exon1", chrom="chr1", start=150, end=250, strand=Strand.POSITIVE))
            self.store.add(GenomicFeature(id="gene2", chrom="chr1", start=300, end=400, strand=Strand.NEGATIVE))
            
            # Chromosome 2
            self.store.add(GenomicFeature(id="gene3", chrom="chr2", start=100, end=200, strand=Strand.POSITIVE))
            self.store.add(GenomicFeature(id="exon2", chrom="chr2", start=250, end=350, strand=Strand.NEGATIVE))
        
        # Create a temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.save_path = Path(self.temp_dir.name) / "feature_store.pkl"

    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()

    def test_save_load(self):
        """Test saving and loading a feature store."""
        # Save the store
        self.store.save(self.save_path)
        
        # Verify the file exists
        self.assertTrue(self.save_path.exists())
        
        # Load the store
        loaded_store = GenomicFeatureStore.load(self.save_path)
        
        # Verify the loaded store has the same type
        self.assertEqual(self.store.store_type, loaded_store.store_type)
        
        # Verify the loaded store has the same chromosomes
        self.assertEqual(self.store.get_chromosomes(), loaded_store.get_chromosomes())
        
        # Verify the loaded store has the same feature count
        self.assertEqual(self.store.get_feature_count(), loaded_store.get_feature_count())
        
        # Verify the loaded store has the same features
        for chrom in self.store.get_chromosomes():
            original_features = self.store.get_by_interval(chrom, 0, 1000)
            loaded_features = loaded_store.get_by_interval(chrom, 0, 1000)
            
            # Check that all features are present
            self.assertEqual(len(original_features), len(loaded_features))
            
            # Check feature details
            for i, (orig, loaded) in enumerate(zip(original_features, loaded_features)):
                self.assertEqual(orig.id, loaded.id)
                self.assertEqual(orig.chrom, loaded.chrom)
                self.assertEqual(orig.start, loaded.start)
                self.assertEqual(orig.end, loaded.end)
                self.assertEqual(orig.strand, loaded.strand)

    def test_save_load_different_store_types(self):
        """Test saving and loading different types of feature stores (except MSI)."""
        # Test all store types except MSI which requires special feature types
        store_types = [StoreType.INTERVAL_TREE, StoreType.BINNED, StoreType.BRUTE_FORCE]
        
        for store_type in store_types:
            # Create a store with the current type
            store = GenomicFeatureStore(store_type=store_type)
            
            # Add a feature
            with store:
                store.add(GenomicFeature(id="test_gene", chrom="chr1", start=100, end=200, strand=Strand.POSITIVE))
            
            # Save the store
            save_path = Path(self.temp_dir.name) / f"feature_store_{store_type.value}.pkl"
            store.save(save_path)
            
            # Load the store
            loaded_store = GenomicFeatureStore.load(save_path)
            
            # Verify the loaded store has the same type
            self.assertEqual(store.store_type, loaded_store.store_type)
            
            # Verify the loaded store has the same feature
            original_features = store.get_by_interval("chr1", 0, 1000)
            loaded_features = loaded_store.get_by_interval("chr1", 0, 1000)
            
            self.assertEqual(len(original_features), len(loaded_features))
            self.assertEqual(original_features[0].chrom, loaded_features[0].chrom)
            self.assertEqual(original_features[0].start, loaded_features[0].start)
            self.assertEqual(original_features[0].end, loaded_features[0].end)
    
    def test_save_load_msi_store(self):
        """Test saving and loading MSI store type."""
        # Create an MSI store
        store = GenomicFeatureStore(store_type=StoreType.MSI)
        
        # Create an MsiSiteRecord
        msi_feature = MsiSiteRecord(
            id="test_msi_site",
            chrom="chr1",
            repeat_unit_length=3,
            repeat_unit_binary=0,  # Placeholder value
            repeat_times=5,
            left_flank_binary=0,   # Placeholder value
            right_flank_binary=0,  # Placeholder value
            repeat_unit_bases="CAG",
            left_flank_bases="ACGT",
            right_flank_bases="TGCA",
            location=100,
            strand=Strand.POSITIVE
        )
        
        # Add the feature to the store
        with store:
            store.add(msi_feature)
        
        # Save the store
        save_path = Path(self.temp_dir.name) / "feature_store_msi.pkl"
        store.save(save_path)
        
        # Load the store
        loaded_store = GenomicFeatureStore.load(save_path)
        
        # Verify the loaded store has the same type
        self.assertEqual(store.store_type, loaded_store.store_type)
        
        # Verify the loaded store has the same feature
        original_features = store.get_by_interval("chr1", 0, 1000)
        loaded_features = loaded_store.get_by_interval("chr1", 0, 1000)
        
        self.assertEqual(len(original_features), len(loaded_features))
        self.assertEqual(original_features[0].chrom, loaded_features[0].chrom)
        self.assertEqual(original_features[0].start, loaded_features[0].start)
        self.assertEqual(original_features[0].end, loaded_features[0].end)

    def test_error_handling(self):
        """Test error handling when loading from non-existent file."""
        non_existent_path = Path(self.temp_dir.name) / "non_existent.pkl"
        
        # Attempt to load from a non-existent file should raise FileNotFoundError
        with self.assertRaises(FileNotFoundError):
            GenomicFeatureStore.load(non_existent_path)


if __name__ == "__main__":
    unittest.main()