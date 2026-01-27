# test_region_chromosome_store.py

import unittest

import numpy as np

from pygnome.feature_store import GenomicFeatureStore
from pygnome.feature_store.region_chromosome_store import (
    RegionFeatureCounter, RegionChromosomeStore, FeatureField
)
from pygnome.genomics import GenomicFeature


class DummyFeature(GenomicFeature):
    def __init__(self, chrom, start, end, label):
        super().__init__(id="", strand=None, chrom=chrom, start=start, end=end)
        self.label = label


# Simple factory for GenomicFeature with label
def region_label_record_factory(chrom, start, end, label):
    return DummyFeature(chrom=chrom, start=start, end=end, label=label)


class TestRegionChromosomeStore(unittest.TestCase):

    def setUp(self):
        self.list_of_test_records = [DummyFeature("chr1", 5000, 5010, 'label1'),
                                     DummyFeature("chr1", 5005, 5010, 'label2'),
                                     DummyFeature("chr2", 600005, 6005010, 'label3'),
                                     ]

    def test_region_feature_counter(self):
        counter = RegionFeatureCounter(bin_size=1000)

        # Add two features in the same bin

        for f in self.list_of_test_records:
            counter.add(f)

        self.assertEqual(counter.get_count("chr1"), 2)
        max_lengths = counter.get_max_lengths("chr1")
        self.assertEqual(list(max_lengths.values()), [10])  # bin 1, max length is 20

    def test_interval_overlap(self):
        """Test finding a feature when the query interval overlaps with it."""
        store, record = self._create_store_with_features()

        features = store.get_by_interval(4590, 5001)
        self.assertEqual(len(features), 1, "Should find feature when interval overlaps with it")
        self.assertEqual(features[0].start, 5000, "Feature start position should match")

    def test_three_records_store_and_interval(self):
        # Counter

        store = self._create_store_with_features()

        # Test get_by_interval
        results = store.get_by_interval('chr1', 4990, 5020)
        self.assertEqual(len(results), 2)
        self.assertEqual(results[0].start, 5000)

        results = store.get_by_interval('chr1', 6000, 7000)
        self.assertEqual(len(results), 0)

    def _create_store_with_features(self):
        counter = RegionFeatureCounter(bin_size=1000)
        for f in self.list_of_test_records:
            counter.add(f)
        store = GenomicFeatureStore()
        count = counter.get_count('chr1')
        max_lengths = counter.get_max_lengths('chr2')
        # define fields for store
        fields = [FeatureField("label", dtype=np.object_)]
        for feat in counter.feature_counts:
            chr_store = RegionChromosomeStore(
                chrom=feat,
                feature_count=count,
                max_lengths_by_bin=max_lengths,
                fields=fields,
                feature_factory=region_label_record_factory,
                bin_size=1000
            )

            chr_store.index_build_start()
            store.chromosomes[feat] = chr_store
        for f in self.list_of_test_records:
            store.add(f)
        for chrom_store in store.chromosomes.values():
            chrom_store.index_build_end()
        return store


if __name__ == "__main__":
    unittest.main()
