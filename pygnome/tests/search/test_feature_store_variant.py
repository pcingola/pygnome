"""Tests for variant search functionality in genomic feature stores."""

import unittest

from pygnome.feature_store import (
    BinnedGenomicStore, GenomicFeatureStore, 
    IntervalTreeStore, StoreType, BruteForceFeatureStore
)
from pygnome.genomics import Strand
from pygnome.genomics.variant import Variant, SNP, Insertion, Deletion


class TestVariantSearch(unittest.TestCase):
    """Test cases for variant search functionality in genomic feature stores."""
    
    def setUp(self):
        """Set up test data."""
        # Create some test variants
        self.variants = []
        
        # Add some SNPs
        for i in range(10):
            pos = i * 100
            self.variants.append(SNP(
                id=f"snp_{i}", 
                chrom="chr1", 
                start=pos, 
                end=pos + 1, 
                strand=Strand.POSITIVE,
                ref="A",
                alt="G"
            ))
        
        # Add some insertions
        for i in range(5):
            pos = i * 200 + 50
            self.variants.append(Insertion(
                id=f"ins_{i}", 
                chrom="chr1", 
                start=pos, 
                end=pos + 1, 
                strand=Strand.POSITIVE,
                ref="A",
                alt="AGT"
            ))
        
        # Add some deletions
        for i in range(5):
            pos = i * 200 + 150
            self.variants.append(Deletion(
                id=f"del_{i}", 
                chrom="chr1", 
                start=pos, 
                end=pos + 3, 
                strand=Strand.POSITIVE,
                ref="ATG",
                alt="A"
            ))
    
    def test_chromosome_store_variant_search(self):
        """Test variant search in chromosome feature stores."""
        for store_class in [BruteForceFeatureStore, IntervalTreeStore, BinnedGenomicStore]:
            # Create store
            store = store_class("chr1")
            
            # Add variants
            with store:
                for variant in self.variants:
                    store.add(variant)
            
            # Test get_by_variant
            # Test SNP
            snp = self.variants[0]
            results = store.get_by_variant(snp.start, snp.ref, snp.alt)
            self.assertEqual(1, len(results), f"{store_class.__name__} failed to find SNP")
            self.assertEqual(snp.id, results[0].id)
            
            # Test insertion
            ins = self.variants[10]
            results = store.get_by_variant(ins.start, ins.ref, ins.alt)
            self.assertEqual(1, len(results), f"{store_class.__name__} failed to find insertion")
            self.assertEqual(ins.id, results[0].id)
            
            # Test deletion
            deletion = self.variants[15]
            results = store.get_by_variant(deletion.start, deletion.ref, deletion.alt)
            self.assertEqual(1, len(results), f"{store_class.__name__} failed to find deletion")
            self.assertEqual(deletion.id, results[0].id)
            
            # Test get_variant
            results = store.get_variant(snp)
            self.assertEqual(1, len(results), f"{store_class.__name__} failed to find variant using get_variant")
            self.assertEqual(snp.id, results[0].id)
    
    def test_genomic_feature_store_variant_search(self):
        """Test variant search in genomic feature stores."""
        for store_type in [StoreType.BRUTE_FORCE, StoreType.INTERVAL_TREE, StoreType.BINNED]:
            # Create store
            store = GenomicFeatureStore(store_type=store_type)
            
            # Add variants
            with store:
                for variant in self.variants:
                    store.add(variant)
            
            # Test get_by_variant
            # Test SNP
            snp = self.variants[0]
            results = store.get_by_variant(snp.chrom, snp.start, snp.ref, snp.alt)
            self.assertEqual(1, len(results), f"{store_type} failed to find SNP")
            self.assertEqual(snp.id, results[0].id)
            
            # Test insertion
            ins = self.variants[10]
            results = store.get_by_variant(ins.chrom, ins.start, ins.ref, ins.alt)
            self.assertEqual(1, len(results), f"{store_type} failed to find insertion")
            self.assertEqual(ins.id, results[0].id)
            
            # Test deletion
            deletion = self.variants[15]
            results = store.get_by_variant(deletion.chrom, deletion.start, deletion.ref, deletion.alt)
            self.assertEqual(1, len(results), f"{store_type} failed to find deletion")
            self.assertEqual(deletion.id, results[0].id)
            
            # Test get_variant
            results = store.get_variant(snp)
            self.assertEqual(1, len(results), f"{store_type} failed to find variant using get_variant")
            self.assertEqual(snp.id, results[0].id)
    
    def test_nonexistent_variant(self):
        """Test searching for a variant that doesn't exist."""
        store = GenomicFeatureStore()
        
        # Add variants
        with store:
            for variant in self.variants:
                store.add(variant)
        
        # Test get_by_variant with nonexistent variant
        results = store.get_by_variant("chr1", 999999, "A", "T")
        self.assertEqual(0, len(results), "Should return empty list for nonexistent variant")
        
        # Test get_variant with nonexistent variant
        nonexistent = SNP(
            id="nonexistent", 
            chrom="chr1", 
            start=999999, 
            end=1000000, 
            strand=Strand.POSITIVE,
            ref="A",
            alt="T"
        )
        results = store.get_variant(nonexistent)
        self.assertEqual(0, len(results), "Should return empty list for nonexistent variant")
        
        # Test get_by_variant with nonexistent chromosome
        results = store.get_by_variant("chr2", 0, "A", "T")
        self.assertEqual(0, len(results), "Should return empty list for nonexistent chromosome")


if __name__ == "__main__":
    unittest.main()