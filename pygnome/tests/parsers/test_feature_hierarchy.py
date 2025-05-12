"""
Tests for the FeatureHierarchy class.
"""

import unittest
from pygnome.parsers.gff.feature_hierarchy import FeatureHierarchy
from pygnome.parsers.gff.record import Record


class TestFeatureHierarchy(unittest.TestCase):
    """Test cases for the FeatureHierarchy class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.hierarchy = FeatureHierarchy()
        
        # Create test records
        self.gene = Record(
            seqid="chr1",
            source="test",
            type="gene",
            start=1000,
            end=5000,
            score=".",
            strand="+",
            phase=".",
            attributes={"ID": "gene1", "Name": "BRCA1"}
        )
        
        self.mrna1 = Record(
            seqid="chr1",
            source="test",
            type="mRNA",
            start=1000,
            end=5000,
            score=".",
            strand="+",
            phase=".",
            attributes={"ID": "mrna1", "Parent": "gene1", "Name": "BRCA1.1"}
        )
        
        self.mrna2 = Record(
            seqid="chr1",
            source="test",
            type="mRNA",
            start=1200,
            end=4800,
            score=".",
            strand="+",
            phase=".",
            attributes={"ID": "mrna2", "Parent": "gene1", "Name": "BRCA1.2"}
        )
        
        self.exon1 = Record(
            seqid="chr1",
            source="test",
            type="exon",
            start=1000,
            end=1200,
            score=".",
            strand="+",
            phase=".",
            attributes={"ID": "exon1", "Parent": "mrna1"}
        )
        
        self.exon2 = Record(
            seqid="chr1",
            source="test",
            type="exon",
            start=1500,
            end=1800,
            score=".",
            strand="+",
            phase=".",
            attributes={"ID": "exon2", "Parent": ["mrna1", "mrna2"]}
        )
        
        self.exon3 = Record(
            seqid="chr1",
            source="test",
            type="exon",
            start=4500,
            end=5000,
            score=".",
            strand="+",
            phase=".",
            attributes={"ID": "exon3", "Parent": "mrna1"}
        )
        
        self.cds1 = Record(
            seqid="chr1",
            source="test",
            type="CDS",
            start=1500,
            end=1800,
            score=".",
            strand="+",
            phase="0",
            attributes={"ID": "cds1", "Parent": "mrna1"}
        )
    
    def test_add_record(self):
        """Test adding records to the hierarchy."""
        # Add records
        self.hierarchy.add_record(self.gene)
        self.hierarchy.add_record(self.mrna1)
        self.hierarchy.add_record(self.mrna2)
        
        # Check ID to record mapping
        self.assertEqual(self.hierarchy.id_to_record["gene1"], self.gene)
        self.assertEqual(self.hierarchy.id_to_record["mrna1"], self.mrna1)
        self.assertEqual(self.hierarchy.id_to_record["mrna2"], self.mrna2)
        
        # Check parent-child relationships
        self.assertEqual(len(self.hierarchy.parent_to_children["gene1"]), 2)
        self.assertIn(self.mrna1, self.hierarchy.parent_to_children["gene1"])
        self.assertIn(self.mrna2, self.hierarchy.parent_to_children["gene1"])
        
        # Check child-parent relationships
        self.assertEqual(len(self.hierarchy.child_to_parents["mrna1"]), 1)
        self.assertEqual(self.hierarchy.child_to_parents["mrna1"][0], "gene1")
    
    def test_build_from_records(self):
        """Test building hierarchy from a list of records."""
        records = [self.gene, self.mrna1, self.mrna2, self.exon1, self.exon2, self.exon3, self.cds1]
        self.hierarchy.build_from_records(records)
        
        # Check parent-child relationships
        self.assertEqual(len(self.hierarchy.parent_to_children["gene1"]), 2)  # 2 mRNAs
        self.assertEqual(len(self.hierarchy.parent_to_children["mrna1"]), 4)  # 3 exons + 1 CDS
        self.assertEqual(len(self.hierarchy.parent_to_children["mrna2"]), 1)  # 1 exon
    
    def test_get_children(self):
        """Test getting children of a feature."""
        records = [self.gene, self.mrna1, self.mrna2, self.exon1, self.exon2, self.exon3, self.cds1]
        self.hierarchy.build_from_records(records)
        
        # Get children of gene
        gene_children = self.hierarchy.get_children("gene1")
        self.assertEqual(len(gene_children), 2)
        self.assertIn(self.mrna1, gene_children)
        self.assertIn(self.mrna2, gene_children)
        
        # Get children of mRNA1
        mrna1_children = self.hierarchy.get_children("mrna1")
        self.assertEqual(len(mrna1_children), 4)
        self.assertIn(self.exon1, mrna1_children)
        self.assertIn(self.exon2, mrna1_children)
        self.assertIn(self.exon3, mrna1_children)
        self.assertIn(self.cds1, mrna1_children)
        
        # Get children of mRNA2
        mrna2_children = self.hierarchy.get_children("mrna2")
        self.assertEqual(len(mrna2_children), 1)
        self.assertIn(self.exon2, mrna2_children)
    
    def test_get_parents(self):
        """Test getting parents of a feature."""
        records = [self.gene, self.mrna1, self.mrna2, self.exon1, self.exon2, self.exon3, self.cds1]
        self.hierarchy.build_from_records(records)
        
        # Get parents of mRNA
        mrna_parents = self.hierarchy.get_parents("mrna1")
        self.assertEqual(len(mrna_parents), 1)
        self.assertEqual(mrna_parents[0], self.gene)
        
        # Get parents of exon with multiple parents
        exon_parents = self.hierarchy.get_parents("exon2")
        self.assertEqual(len(exon_parents), 2)
        self.assertIn(self.mrna1, exon_parents)
        self.assertIn(self.mrna2, exon_parents)
    
    def test_get_descendants(self):
        """Test getting all descendants of a feature."""
        records = [self.gene, self.mrna1, self.mrna2, self.exon1, self.exon2, self.exon3, self.cds1]
        self.hierarchy.build_from_records(records)
        
        # Get all descendants of gene
        descendants = self.hierarchy.get_descendants("gene1")
        self.assertEqual(len(descendants), 6)  # 2 mRNAs + 3 exons + 1 CDS
        
        # Test with max_depth
        descendants = self.hierarchy.get_descendants("gene1", max_depth=1)
        self.assertEqual(len(descendants), 2)  # Only direct children (mRNAs)
    
    def test_get_ancestors(self):
        """Test getting all ancestors of a feature."""
        records = [self.gene, self.mrna1, self.mrna2, self.exon1, self.exon2, self.exon3, self.cds1]
        self.hierarchy.build_from_records(records)
        
        # Get ancestors of exon
        ancestors = self.hierarchy.get_ancestors("exon1")
        self.assertEqual(len(ancestors), 2)  # mRNA1 + gene
        
        # Test with max_depth
        ancestors = self.hierarchy.get_ancestors("exon1", max_depth=1)
        self.assertEqual(len(ancestors), 1)  # Only direct parent (mRNA1)


if __name__ == "__main__":
    unittest.main()