"""
Tests for the GFF3 parser.
"""

import unittest
from pathlib import Path
from pygnome.parsers.gff3_parser import Gff3Parser
from pygnome.genomics.strand import Strand


class TestGff3Parser(unittest.TestCase):
    """Test cases for the GFF3 parser."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.parser = Gff3Parser()
        self.test_file = Path(__file__).parent.parent / "data" / "test.gff3"
    
    def test_parse_attributes(self):
        """Test parsing GFF3 attributes."""
        # Simple attributes
        attr_string = "ID=gene00001;Name=EDEN"
        attributes = self.parser._parse_attributes(attr_string)
        self.assertEqual(attributes["ID"], "gene00001")
        self.assertEqual(attributes["Name"], "EDEN")
        
        # Multiple values for an attribute
        attr_string = "ID=exon00002;Parent=mRNA00001,mRNA00002"
        attributes = self.parser._parse_attributes(attr_string)
        self.assertEqual(attributes["ID"], "exon00002")
        self.assertEqual(attributes["Parent"], ["mRNA00001", "mRNA00002"])
        
        # URL-encoded values
        attr_string = "ID=gene%3A001;Note=Complex%20note%3B%20with%20semicolons"
        attributes = self.parser._parse_attributes(attr_string)
        self.assertEqual(attributes["ID"], "gene:001")
        self.assertEqual(attributes["Note"], "Complex note; with semicolons")
        
        # Empty attribute string
        attributes = self.parser._parse_attributes(".")
        self.assertEqual(attributes, {})
    
    def test_parse_line(self):
        """Test parsing a GFF3 line."""
        line = "ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN"
        record = self.parser._parse_line(line)
        
        self.assertEqual(record.seqid, "ctg123")
        self.assertEqual(record.source, ".")
        self.assertEqual(record.type, "gene")
        self.assertEqual(record.start, 1000)
        self.assertEqual(record.end, 9000)
        self.assertIsNone(record.score)
        self.assertEqual(record.strand, Strand.POSITIVE)
        self.assertIsNone(record.phase)
        self.assertEqual(record.attributes["ID"], "gene00001")
        self.assertEqual(record.attributes["Name"], "EDEN")
        
        # Invalid line (too few fields)
        line = "ctg123\t.\tgene\t1000\t9000\t.\t+"
        record = self.parser._parse_line(line)
        self.assertIsNone(record)
    
    def test_parse_file(self):
        """Test parsing a GFF3 file."""
        records = list(self.parser.parse(self.test_file))
        
        # Check number of records
        self.assertEqual(len(records), 14)
        
        # Check first record
        first_record = records[0]
        self.assertEqual(first_record.seqid, "ctg123")
        self.assertEqual(first_record.type, "gene")
        self.assertEqual(first_record.start, 1000)
        self.assertEqual(first_record.end, 9000)
        self.assertEqual(first_record.attributes["ID"], "gene00001")
        self.assertEqual(first_record.attributes["Name"], "EDEN")
        
        # Check a record with multiple parents
        exon_record = None
        for record in records:
            if record.type == "exon" and record.get_attribute("ID") == "exon00004":
                exon_record = record
                break
        
        self.assertIsNotNone(exon_record)
        self.assertEqual(exon_record.start, 5000)
        self.assertEqual(exon_record.end, 5500)
        self.assertEqual(exon_record.get_attribute("Parent"), 
                         ["mRNA00001", "mRNA00002", "mRNA00003"])


if __name__ == "__main__":
    unittest.main()