"""
Tests for the VcfFormat class.
"""
import unittest

from pygnome.parsers.vcf.vcf_header import VcfHeader
from pygnome.parsers.vcf.vcf_format import VcfFormat


class TestVcfFormat(unittest.TestCase):
    """Test cases for the VcfFormat class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a header
        self.header = VcfHeader()
        
        # Add meta lines
        self.header.add_meta_line('##fileformat=VCFv4.5')
        self.header.add_meta_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        self.header.add_meta_line('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">')
        self.header.add_meta_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">')
        self.header.add_meta_line('##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">')
        
        # Create a format string
        self.format_str = "GT:GQ:DP:HQ"
        self.format = VcfFormat(self.format_str, self.header)
    
    def test_get_keys(self):
        """Test getting the keys in the FORMAT field."""
        keys = self.format.get_keys()
        self.assertEqual(keys, ["GT", "GQ", "DP", "HQ"])
    
    def test_has_key(self):
        """Test checking if a key is present in the FORMAT field."""
        self.assertTrue(self.format.has_key("GT"))
        self.assertTrue(self.format.has_key("GQ"))
        self.assertTrue(self.format.has_key("DP"))
        self.assertTrue(self.format.has_key("HQ"))
        self.assertFalse(self.format.has_key("XX"))
    
    def test_add_key(self):
        """Test adding a key to the FORMAT field."""
        # Add a new key
        self.header.add_meta_line('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">')
        self.format.add_key("PL")
        
        # Check that the key was added
        self.assertTrue(self.format.has_key("PL"))
        self.assertEqual(self.format.get_keys(), ["GT", "GQ", "DP", "HQ", "PL"])
        
        # Check the string representation
        self.assertEqual(str(self.format), "GT:GQ:DP:HQ:PL")
    
    def test_remove_key(self):
        """Test removing a key from the FORMAT field."""
        # Remove a key
        self.format.remove_key("HQ")
        
        # Check that the key was removed
        self.assertFalse(self.format.has_key("HQ"))
        self.assertEqual(self.format.get_keys(), ["GT", "GQ", "DP"])
        
        # Check the string representation
        self.assertEqual(str(self.format), "GT:GQ:DP")
    
    def test_to_string(self):
        """Test converting the FORMAT field to a string."""
        self.assertEqual(self.format.to_string(), "GT:GQ:DP:HQ")
        
        # Test with an empty format
        empty_format = VcfFormat("", self.header)
        self.assertEqual(empty_format.to_string(), ".")
    
    def test_empty_format(self):
        """Test handling of empty FORMAT fields."""
        empty_format = VcfFormat(".", self.header)
        self.assertEqual(empty_format.get_keys(), [])
        self.assertFalse(empty_format.has_key("GT"))
        self.assertEqual(empty_format.to_string(), ".")


if __name__ == '__main__':
    unittest.main()