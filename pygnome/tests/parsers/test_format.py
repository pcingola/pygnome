"""
Tests for the Format enum.
"""

import unittest
from pygnome.parsers.gff.format import Format


class TestFormat(unittest.TestCase):
    """Test cases for the Format enum."""
    
    def test_format_values(self):
        """Test Format enum values."""
        self.assertEqual(Format.GFF2.value, "gff2")
        self.assertEqual(Format.GFF3.value, "gff3")
        self.assertEqual(Format.GTF.value, "gtf")
    
    def test_format_from_string(self):
        """Test creating Format from string values."""
        self.assertEqual(Format("gff2"), Format.GFF2)
        self.assertEqual(Format("gff3"), Format.GFF3)
        self.assertEqual(Format("gtf"), Format.GTF)
        
        # Test with invalid value
        with self.assertRaises(ValueError):
            Format("invalid")
    
    def test_format_string_representation(self):
        """Test string representation of Format."""
        self.assertEqual(str(Format.GFF2), "gff2")
        self.assertEqual(str(Format.GFF3), "gff3")
        self.assertEqual(str(Format.GTF), "gtf")


if __name__ == "__main__":
    unittest.main()