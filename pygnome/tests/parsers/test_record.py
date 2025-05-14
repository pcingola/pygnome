
"""
Tests for the Record class.
"""

import unittest
from pygnome.parsers.gff.record import Record
from pygnome.genomics.strand import Strand


class TestRecord(unittest.TestCase):
    """Test cases for the Record class."""
    
    def test_record_initialization(self):
        """Test basic Record initialization."""
        record = Record(
            seqid="chr1",
            source="test",
            type="gene",
            start="1000",
            end="2000",
            score=".",
            strand="+",
            phase=".",
            attributes={"ID": "gene1", "Name": "BRCA1"}
        )
        
        # Check field conversions
        self.assertEqual(record.seqid, "chr1")
        self.assertEqual(record.source, "test")
        self.assertEqual(record.type, "gene")
        self.assertEqual(record.start, 1000)  # Converted to int
        self.assertEqual(record.end, 2000)    # Converted to int
        self.assertIsNone(record.score)       # "." converted to None
        self.assertEqual(record.strand, Strand.POSITIVE)  # "+" converted to Strand.POSITIVE
        self.assertIsNone(record.phase)       # "." converted to None
        self.assertEqual(record.attributes["ID"], "gene1")
        self.assertEqual(record.attributes["Name"], "BRCA1")
    
    def test_record_with_numeric_values(self):
        """Test Record with numeric values."""
        record = Record(
            seqid="chr1",
            source="test",
            type="gene",
            start=1000,
            end=2000,
            score=95.6,
            strand="-",
            phase=0,
            attributes={}
        )
        
        self.assertEqual(record.start, 1000)
        self.assertEqual(record.end, 2000)
        self.assertEqual(record.score, 95.6)
        self.assertEqual(record.strand, Strand.NEGATIVE)
        self.assertEqual(record.phase, 0)
    
    def test_get_set_attribute(self):
        """Test get_attribute and set_attribute methods."""
        record = Record(
            seqid="chr1",
            source="test",
            type="gene",
            start=1000,
            end=2000,
            score=".",
            strand=".",
            phase=".",
            attributes={"ID": "gene1"}
        )
        
        # Test get_attribute
        self.assertEqual(record.get_attribute("ID"), "gene1")
        self.assertIsNone(record.get_attribute("Name"))
        self.assertEqual(record.get_attribute("Name", "unknown"), "unknown")
        
        # Test set_attribute
        record.set_attribute("Name", "BRCA1")
        self.assertEqual(record.get_attribute("Name"), "BRCA1")
    
    def test_string_representation(self):
        """Test string representation of Record."""
        record = Record(
            seqid="chr1",
            source="test",
            type="gene",
            start=1000,
            end=2000,
            score=95.6,
            strand="+",
            phase=0,
            attributes={"ID": "gene1", "Name": "BRCA1"}
        )
        
        # Convert to string and check format
        record_str = str(record)
        self.assertIn("chr1", record_str)
        self.assertIn("test", record_str)
        self.assertIn("gene", record_str)
        self.assertIn("1000", record_str)
        self.assertIn("2000", record_str)
        self.assertIn("95.6", record_str)
        self.assertIn("+", record_str)
        self.assertIn("0", record_str)
        self.assertTrue("ID=gene1" in record_str)
        self.assertTrue("Name=BRCA1" in record_str)


if __name__ == "__main__":
    unittest.main()