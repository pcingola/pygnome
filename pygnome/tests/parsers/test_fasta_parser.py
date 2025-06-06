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
        # Use absolute path to ensure file can be found regardless of working directory
        current_dir = Path(__file__).parent
        self.test_file = current_dir.parent / "data" / "sequences" / "test.fasta"
    
    def test_parse(self):
        """Test parsing a FASTA file."""
        records = FastaParser(self.test_file).load()
        
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
        sequences = FastaParser(self.test_file).load_as_dict()
        
        self.assertEqual(len(sequences), 3)
        self.assertEqual(sequences["seq1"].sequence, "ACGTACGTACGTACGT")
        self.assertEqual(sequences["seq2"].sequence, "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT")
        self.assertEqual(sequences["seq3"].sequence, "ACGTNRYSWKMBDHV")

    def test_parse_as_dna_strings(self):
        """Test parsing as DnaString objects."""
        sequences = FastaParser(self.test_file, length_convert_to_dna_string=-1).load_as_dict()
        print(sequences)

        self.assertEqual(len(sequences), 3)
        self.assertIsInstance(sequences["seq1"].sequence, DnaString)
        self.assertEqual(str(sequences["seq1"].sequence), "ACGTACGTACGTACGT")
    
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