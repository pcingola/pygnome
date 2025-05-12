"""
Tests for DnaString and RnaString classes.
"""

import unittest

from pygnome.sequences.dna_string import DnaString
from pygnome.sequences.rna_string import RnaString


class TestDnaString(unittest.TestCase):
    """Test cases for the DnaString class."""
    
    def test_initialization(self):
        """Test initialization of DnaString."""
        dna = DnaString("ACGT")
        self.assertEqual(len(dna), 4)
        self.assertEqual(str(dna), "ACGT")
    
    def test_getitem(self):
        """Test accessing individual nucleotides."""
        dna = DnaString("ACGT")
        self.assertEqual(dna[0], "A")
        self.assertEqual(dna[1], "C")
        self.assertEqual(dna[2], "G")
        self.assertEqual(dna[3], "T")
        
        # Test negative indices
        self.assertEqual(dna[-1], "T")
        self.assertEqual(dna[-2], "G")
    
    def test_slicing(self):
        """Test slicing operations."""
        dna = DnaString("ACGTACGT")
        self.assertEqual(dna[1:3], "CG")
        self.assertEqual(dna[::2], "AGAG")
        self.assertEqual(dna[::-1], "TGCATGCA")
    
    def test_substring(self):
        """Test substring method."""
        dna = DnaString("ACGTACGT")
        self.assertEqual(dna.substring(2, 4), "GTAC")
        self.assertEqual(dna.substring(0, 4), "ACGT")
    
    def test_equality(self):
        """Test equality comparison."""
        dna1 = DnaString("ACGT")
        dna2 = DnaString("ACGT")
        dna3 = DnaString("ACGTA")
        
        self.assertEqual(dna1, dna2)
        self.assertNotEqual(dna1, dna3)
    
    def test_ambiguous_nucleotides(self):
        """Test handling of ambiguous nucleotides."""
        dna = DnaString("ACGTN")
        # N should be converted to A in the 2-bit representation
        self.assertEqual(dna[4], "A")


class TestRnaString(unittest.TestCase):
    """Test cases for the RnaString class."""
    
    def test_initialization(self):
        """Test initialization of RnaString."""
        rna = RnaString("ACGU")
        self.assertEqual(len(rna), 4)
        self.assertEqual(str(rna), "ACGU")
    
    def test_t_to_u_conversion(self):
        """Test automatic conversion of T to U."""
        rna = RnaString("ACGT")
        self.assertEqual(str(rna), "ACGU")
    
    def test_getitem(self):
        """Test accessing individual nucleotides."""
        rna = RnaString("ACGU")
        self.assertEqual(rna[0], "A")
        self.assertEqual(rna[1], "C")
        self.assertEqual(rna[2], "G")
        self.assertEqual(rna[3], "U")
        
        # Test negative indices
        self.assertEqual(rna[-1], "U")
        self.assertEqual(rna[-2], "G")
    
    def test_slicing(self):
        """Test slicing operations."""
        rna = RnaString("ACGUACGU")
        self.assertEqual(rna[1:3], "CG")
        self.assertEqual(rna[::2], "AGAG")
        self.assertEqual(rna[::-1], "UGCAUGCA")
    
    def test_substring(self):
        """Test substring method."""
        rna = RnaString("ACGUACGU")
        self.assertEqual(rna.substring(2, 4), "GUAC")
        self.assertEqual(rna.substring(0, 4), "ACGU")
    
    def test_equality(self):
        """Test equality comparison."""
        rna1 = RnaString("ACGU")
        rna2 = RnaString("ACGU")
        rna3 = RnaString("ACGUA")
        
        self.assertEqual(rna1, rna2)
        self.assertNotEqual(rna1, rna3)


if __name__ == "__main__":
    unittest.main()