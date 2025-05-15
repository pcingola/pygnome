"""
Tests for the DnaStringArray class.
"""

import unittest
import numpy as np

from pygnome.sequences.dna_string import DnaString
from pygnome.sequences.dna_string_array import DnaStringArray


class TestDnaStringArray(unittest.TestCase):
    """Test cases for the DnaStringArray class."""
    
    def test_init(self):
        """Test initialization with different capacities."""
        # Default initialization
        array = DnaStringArray()
        self.assertEqual(len(array), 0)
        
        # Custom initialization
        array = DnaStringArray(initial_data_bytes=8000, initial_strings=2000)
        self.assertEqual(len(array), 0)
    
    def test_add_single(self):
        """Test adding a single sequence."""
        array = DnaStringArray()
        
        # Add a sequence
        idx = array.add("ACGT")
        self.assertEqual(idx, 0)
        self.assertEqual(len(array), 1)
        self.assertEqual(array[0], "ACGT")
        
        # Add another sequence
        idx = array.add("TGCA")
        self.assertEqual(idx, 1)
        self.assertEqual(len(array), 2)
        self.assertEqual(array[0], "ACGT")
        self.assertEqual(array[1], "TGCA")
    
    def test_add_multiple(self):
        """Test adding multiple sequences at once."""
        array = DnaStringArray()
        
        # Add multiple sequences
        sequences = ["ACGT", "TGCA", "AAAA", "CCCC", "GGGG", "TTTT"]
        indices = array.add_multiple(sequences)
        
        self.assertEqual(indices, [0, 1, 2, 3, 4, 5])
        self.assertEqual(len(array), 6)
        
        # Verify all sequences
        for i, seq in enumerate(sequences):
            self.assertEqual(array[i], seq)
    
    def test_empty_sequence(self):
        """Test handling of empty sequences."""
        array = DnaStringArray()
        
        # Add an empty sequence
        idx = array.add("")
        self.assertEqual(idx, 0)
        self.assertEqual(len(array), 1)
        self.assertEqual(array[0], "")
        self.assertEqual(array.get_length(0), 0)
        
        # Add a non-empty sequence after an empty one
        idx = array.add("ACGT")
        self.assertEqual(idx, 1)
        self.assertEqual(len(array), 2)
        self.assertEqual(array[0], "")
        self.assertEqual(array[1], "ACGT")
    
    def test_long_sequence(self):
        """Test handling of sequences longer than 16 nucleotides."""
        array = DnaStringArray()
        
        # Add a sequence that spans multiple integers
        long_seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"  # 32 nucleotides
        idx = array.add(long_seq)
        
        self.assertEqual(idx, 0)
        self.assertEqual(len(array), 1)
        self.assertEqual(array[0], long_seq)
        self.assertEqual(array.get_length(0), len(long_seq))
    
    def test_get_subsequence(self):
        """Test extracting subsequences."""
        array = DnaStringArray()
        
        # Add a sequence
        seq = "ACGTACGTACGTACGT"  # 16 nucleotides
        idx = array.add(seq)
        
        # Extract various subsequences
        self.assertEqual(array.get_subsequence(0, 0, 4), "ACGT")
        self.assertEqual(array.get_subsequence(0, 4, 4), "ACGT")
        self.assertEqual(array.get_subsequence(0, 8, 4), "ACGT")
        self.assertEqual(array.get_subsequence(0, 12, 4), "ACGT")
        
        # Extract with default length (to end)
        self.assertEqual(array.get_subsequence(0, 12), "ACGT")
        
        # Extract with negative start
        self.assertEqual(array.get_subsequence(0, -4), "ACGT")
    
    def test_to_dna_string(self):
        """Test conversion to DnaString objects."""
        array = DnaStringArray()
        
        # Add a sequence
        seq = "ACGTACGT"
        idx = array.add(seq)
        
        # Convert to DnaString
        dna_string = array.to_dna_string(0)
        
        self.assertIsInstance(dna_string, DnaString)
        self.assertEqual(str(dna_string), seq)
    
    def test_get_stats(self):
        """Test getting statistics about memory usage."""
        array = DnaStringArray()
        
        # Empty array
        count, total_nt, bits_per_nt = array.get_stats()
        self.assertEqual(count, 0)
        self.assertEqual(total_nt, 0)
        self.assertEqual(bits_per_nt, 0)
        
        # Add some sequences
        array.add("ACGT")  # 4 nucleotides
        array.add("TGCA")  # 4 nucleotides
        
        count, total_nt, bits_per_nt = array.get_stats()
        self.assertEqual(count, 2)
        self.assertEqual(total_nt, 8)
        # bits_per_nt will be slightly more than 2 due to overhead
        self.assertGreater(bits_per_nt, 2.0)
    
    def test_large_number_of_sequences(self):
        """Test handling a large number of sequences."""
        # Create with larger initial capacity
        array = DnaStringArray(initial_data_bytes=10000, initial_strings=1000)
        
        # Add 1000 short sequences
        sequences = ["ACGT" for _ in range(1000)]
        indices = array.add_multiple(sequences)
        
        self.assertEqual(len(indices), 1000)
        self.assertEqual(len(array), 1000)
        
        # Verify some random sequences
        for i in range(0, 1000, 100):
            self.assertEqual(array[i], "ACGT")
    
    def test_error_handling(self):
        """Test error handling for invalid operations."""
        array = DnaStringArray()
        
        # Add a sequence
        array.add("ACGT")
        
        # Test index out of range
        with self.assertRaises(IndexError):
            array[-2]  # Negative index out of range
        
        with self.assertRaises(IndexError):
            array[1]  # Positive index out of range
        
        # Test invalid subsequence parameters
        with self.assertRaises(IndexError):
            array.get_subsequence(0, 5, 1)  # Start out of range
        
        with self.assertRaises(IndexError):
            array.get_subsequence(0, 0, 5)  # Length out of range


if __name__ == "__main__":
    unittest.main()