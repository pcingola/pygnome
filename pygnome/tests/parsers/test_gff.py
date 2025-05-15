"""
Tests for the Gff context manager.
"""

import unittest
from pathlib import Path
from pygnome.parsers.gff import Gff
from pygnome.parsers.gff.format import Format
from pygnome.genomics.strand import Strand


class TestGff(unittest.TestCase):
    """Test cases for the Gff context manager."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(__file__).parent.parent / "data"
        self.gff3_file = self.test_dir / "test.gff3"
        self.gtf_file = self.test_dir / "test.gtf"
        self.gff2_file = self.test_dir / "test.gff"
    
    def test_format_detection(self):
        """Test format auto-detection."""
        with Gff(self.gff3_file) as gff:
            self.assertEqual(gff.format, Format.GFF3)
        
        with Gff(self.gtf_file) as gff:
            self.assertEqual(gff.format, Format.GTF)
        
        with Gff(self.gff2_file) as gff:
            self.assertEqual(gff.format, Format.GFF2)
    
    def test_explicit_format(self):
        """Test explicitly specifying format."""
        with Gff(self.gff3_file, format=Format.GFF3) as gff:
            self.assertEqual(gff.format, Format.GFF3)
        
        # Test with incorrect format (should still work but parse differently)
        with Gff(self.gff3_file, format=Format.GTF) as gff:
            self.assertEqual(gff.format, Format.GTF)
    
    def test_iteration(self):
        """Test iterating over records."""
        with Gff(self.gff3_file) as gff:
            records = list(gff)
            self.assertEqual(len(records), 14)
            
            # Check first record
            first = records[0]
            self.assertEqual(first.chrom, "ctg123")
            self.assertEqual(first.type, "gene")
            self.assertEqual(first.start, 1000)
            self.assertEqual(first.end, 9000)
    
    def test_get_features_by_type(self):
        """Test filtering features by type."""
        with Gff(self.gff3_file) as gff:
            # Consume iterator to populate records
            list(gff)
            
            # Get features by type
            genes = gff.get_features_by_type("gene")
            mrnas = gff.get_features_by_type("mRNA")
            exons = gff.get_features_by_type("exon")
            
            self.assertEqual(len(genes), 1)
            self.assertEqual(len(mrnas), 3)
            self.assertEqual(len(exons), 5)
    
    def test_get_features_by_location(self):
        """Test filtering features by location."""
        with Gff(self.gff3_file) as gff:
            # Consume iterator to populate records
            list(gff)
            
            # Get features in a specific region
            features = gff.get_features_by_location("ctg123", 1000, 1500)
            
            # Should include gene, TF_binding_site, mRNAs, exons, and CDS
            self.assertGreaterEqual(len(features), 7)
            
            # Check types
            types = set(f.type for f in features)
            self.assertIn("gene", types)
            self.assertIn("mRNA", types)
            self.assertIn("exon", types)
    
    def test_get_features_by_attribute(self):
        """Test filtering features by attribute."""
        with Gff(self.gff3_file) as gff:
            # Consume iterator to populate records
            list(gff)
            
            # Get features with specific ID
            features = gff.get_features_by_attribute("ID", "gene00001")
            self.assertEqual(len(features), 1)
            self.assertEqual(features[0].type, "gene")
            
            # Get features with any Parent attribute
            features = gff.get_features_by_attribute("Parent")
            self.assertGreaterEqual(len(features), 13)  # All except gene
            
            # Get features with specific Parent
            features = gff.get_features_by_attribute("Parent", "gene00001")
            self.assertEqual(len(features), 4)  # TF_binding_site + 3 mRNAs
    
    def test_build_hierarchy(self):
        """Test building feature hierarchy."""
        with Gff(self.gff3_file) as gff:
            # Consume iterator to populate records
            list(gff)
            
            # Build hierarchy
            hierarchy = gff.build_hierarchy()
            
            # Get children of gene
            gene_children = hierarchy.get_children("gene00001")
            self.assertEqual(len(gene_children), 4)  # TF_binding_site + 3 mRNAs
            
            # Get children of mRNA
            mrna_children = hierarchy.get_children("mRNA00001")
            self.assertGreaterEqual(len(mrna_children), 4)  # exons + CDS
            
            # Get parents of exon
            exon_parents = hierarchy.get_parents("exon00004")
            self.assertEqual(len(exon_parents), 3)  # 3 mRNAs


if __name__ == "__main__":
    unittest.main()