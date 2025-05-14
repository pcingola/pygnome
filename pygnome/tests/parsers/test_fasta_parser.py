"""
Tests for FASTA parser.
"""

import unittest
from pathlib import Path

from pygnome.parsers.fasta.fasta_parser import FastaParser, FastaRecord
from pygnome.sequences.dna_string import DnaString
from pygnome.sequences.rna_string import RnaString


class TestFastaParser(unittest.TestCase):
    """Test cases for the FASTA parser."""
    
    def setUp(self):
        """Set up test data."""
        self.test_file = Path("pygnome/tests/data/sequences/test.fasta")
    
    def test_parse(self):
        """Test parsing a FASTA file."""
        records = list(FastaParser.parse(self.test_file))
        
        self.assertEqual(len(records), 3)
        
        # Check first record
        self.assertEqual(records[0].identifier, "seq1")
        self.assertEqual(records[0].description, "Test sequence 1")
        self.assertEqual(records[0].sequence, "ACGTACGTACGTACGT")
        
        # Check second record
        self.assertEqual(records[1].identifier, "seq2")
        self.assertEqual(records[1].description, "Test sequence 2")
        self.assertEqual(records[1].sequence, "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT")
        
        # Check third record
        self.assertEqual(records[2].identifier, "seq3")
        self.assertEqual(records[2].description, "Test sequence with ambiguous nucleotides")
        self.assertEqual(records[2].sequence, "ACGTNRYSWKMBDHV")
    
    def test_context_manager(self):
        """Test using the parser as a context manager."""
        with FastaParser(self.test_file) as parser:
            records = list(parser)
            
            self.assertEqual(len(records), 3)
            self.assertEqual(records[0].identifier, "seq1")
    
    def test_parse_as_dict(self):
        """Test parsing a FASTA file as a dictionary."""
        sequences = FastaParser.parse_as_dict(self.test_file)
        
        self.assertEqual(len(sequences), 3)
        self.assertEqual(sequences["seq1"], "ACGTACGTACGTACGT")
        self.assertEqual(sequences["seq2"], "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT")
        self.assertEqual(sequences["seq3"], "ACGTNRYSWKMBDHV")
    
    def test_parse_first(self):
        """Test parsing only the first record."""
        record = FastaParser.parse_first(self.test_file)
        
        self.assertIsNotNone(record)
        self.assertEqual(record.identifier, "seq1")
        self.assertEqual(record.description, "Test sequence 1")
        self.assertEqual(record.sequence, "ACGTACGTACGTACGT")
    
    def test_parse_as_dna_strings(self):
        """Test parsing as DnaString objects."""
        sequences = FastaParser.parse_as_dna_strings(self.test_file)
        
        self.assertEqual(len(sequences), 3)
        self.assertIsInstance(sequences["seq1"], DnaString)
        self.assertEqual(str(sequences["seq1"]), "ACGTACGTACGTACGT")
    
    def test_parse_as_rna_strings(self):
        """Test parsing as RnaString objects."""
        sequences = FastaParser.parse_as_rna_strings(self.test_file)
        
        self.assertEqual(len(sequences), 3)
        self.assertIsInstance(sequences["seq1"], RnaString)
        self.assertEqual(str(sequences["seq1"]), "ACGUACGUACGUACGU")  # Note: T converted to U
    
    def test_fasta_record_str(self):
        """Test string representation of FastaRecord."""
        record = FastaRecord("test_id", "ACGT", "Test description")
        expected = ">test_id Test description\nACGT"
        self.assertEqual(str(record), expected)
        
        # Test long sequence formatting (80 chars per line)
        long_seq = "A" * 100
        record = FastaRecord("long", long_seq)
        expected = ">long\n" + "A" * 80 + "\n" + "A" * 20
        self.assertEqual(str(record), expected)


if __name__ == "__main__":
    unittest.main()