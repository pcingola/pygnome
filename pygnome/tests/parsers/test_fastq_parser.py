"""
Tests for FASTQ parser.
"""

import unittest
from pathlib import Path

from pygnome.parsers.fastq_parser import FastqParser, FastqRecord
from pygnome.sequences.dna_string import DnaString
from pygnome.sequences.rna_string import RnaString


class TestFastqParser(unittest.TestCase):
    """Test cases for the FASTQ parser."""
    
    def setUp(self):
        """Set up test data."""
        self.test_file = Path("tests/data/sequences/test.fastq")
    
    def test_parse(self):
        """Test parsing a FASTQ file."""
        records = list(FastqParser.parse(self.test_file))
        
        self.assertEqual(len(records), 3)
        
        # Check first record
        self.assertEqual(records[0].identifier, "read1")
        self.assertEqual(records[0].description, "Test read 1")
        self.assertEqual(records[0].sequence, "ACGTACGTACGTACGT")
        self.assertEqual(records[0].quality, "IIIIIIIIIIIIIIII")
        
        # Check second record
        self.assertEqual(records[1].identifier, "read2")
        self.assertEqual(records[1].description, "Test read 2")
        self.assertEqual(records[1].sequence, "AAAACCCCGGGGTTTT")
        self.assertEqual(records[1].quality, "HHHHIIIIJJJJKKKK")
        
        # Check third record
        self.assertEqual(records[2].identifier, "read3")
        self.assertEqual(records[2].description, "Test read with ambiguous nucleotides")
        self.assertEqual(records[2].sequence, "ACGTNRYSWKMBDHV")
        self.assertEqual(records[2].quality, "IIIIHHHHGGGGFFFF")
    
    def test_context_manager(self):
        """Test using the parser as a context manager."""
        with FastqParser(self.test_file) as parser:
            records = list(parser)
            
            self.assertEqual(len(records), 3)
            self.assertEqual(records[0].identifier, "read1")
    
    def test_parse_as_dict(self):
        """Test parsing a FASTQ file as a dictionary."""
        sequences = FastqParser.parse_as_dict(self.test_file)
        
        self.assertEqual(len(sequences), 3)
        self.assertEqual(sequences["read1"], ("ACGTACGTACGTACGT", "IIIIIIIIIIIIIIII"))
        self.assertEqual(sequences["read2"], ("AAAACCCCGGGGTTTT", "HHHHIIIIJJJJKKKK"))
        self.assertEqual(sequences["read3"], ("ACGTNRYSWKMBDHV", "IIIIHHHHGGGGFFFF"))
    
    def test_parse_first(self):
        """Test parsing only the first record."""
        record = FastqParser.parse_first(self.test_file)
        
        self.assertIsNotNone(record)
        self.assertEqual(record.identifier, "read1")
        self.assertEqual(record.description, "Test read 1")
        self.assertEqual(record.sequence, "ACGTACGTACGTACGT")
        self.assertEqual(record.quality, "IIIIIIIIIIIIIIII")
    
    def test_parse_as_dna_strings(self):
        """Test parsing as DnaString objects."""
        sequences = FastqParser.parse_as_dna_strings(self.test_file)
        
        self.assertEqual(len(sequences), 3)
        dna_string, quality = sequences["read1"]
        self.assertIsInstance(dna_string, DnaString)
        self.assertEqual(str(dna_string), "ACGTACGTACGTACGT")
        self.assertEqual(quality, "IIIIIIIIIIIIIIII")
    
    def test_parse_as_rna_strings(self):
        """Test parsing as RnaString objects."""
        sequences = FastqParser.parse_as_rna_strings(self.test_file)
        
        self.assertEqual(len(sequences), 3)
        rna_string, quality = sequences["read1"]
        self.assertIsInstance(rna_string, RnaString)
        self.assertEqual(str(rna_string), "ACGUACGUACGUACGU")  # Note: T converted to U
        self.assertEqual(quality, "IIIIIIIIIIIIIIII")
    
    def test_fastq_record_str(self):
        """Test string representation of FastqRecord."""
        record = FastqRecord("test_id", "ACGT", "IIII", "Test description")
        expected = "@test_id Test description\nACGT\n+\nIIII"
        self.assertEqual(str(record), expected)
    
    def test_get_quality_scores(self):
        """Test converting quality strings to Phred scores."""
        record = FastqRecord("test_id", "ACGT", "IIII")
        # 'I' has ASCII value 73, so Phred score is 73 - 33 = 40
        self.assertEqual(record.get_quality_scores(), [40, 40, 40, 40])
    
    def test_validation(self):
        """Test validation of sequence and quality length."""
        # Should not raise an error
        FastqRecord("test_id", "ACGT", "IIII")
        
        # Should raise an error
        with self.assertRaises(ValueError):
            FastqRecord("test_id", "ACGT", "III")


if __name__ == "__main__":
    unittest.main()