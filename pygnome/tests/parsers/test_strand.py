"""
Tests for the Strand enum.
"""

import unittest
from pygnome.genomics.strand import Strand


class TestStrand(unittest.TestCase):
    """Test cases for the Strand enum."""
    
    def test_strand_values(self):
        """Test Strand enum values."""
        self.assertEqual(Strand.POSITIVE.value, "+")
        self.assertEqual(Strand.NEGATIVE.value, "-")
        self.assertEqual(Strand.UNSTRANDED.value, ".")
        self.assertEqual(Strand.UNKNOWN.value, "?")
    
    def test_strand_from_string(self):
        """Test creating Strand from string values."""
        self.assertEqual(Strand("+"), Strand.POSITIVE)
        self.assertEqual(Strand("-"), Strand.NEGATIVE)
        self.assertEqual(Strand("."), Strand.UNSTRANDED)
        self.assertEqual(Strand("?"), Strand.UNKNOWN)
        
        # Test with invalid value
        with self.assertRaises(ValueError):
            Strand("invalid")
    
    def test_strand_string_representation(self):
        """Test string representation of Strand."""
        self.assertEqual(str(Strand.POSITIVE), "+")
        self.assertEqual(str(Strand.NEGATIVE), "-")
        self.assertEqual(str(Strand.UNSTRANDED), ".")
        self.assertEqual(str(Strand.UNKNOWN), "?")


if __name__ == "__main__":
    unittest.main()