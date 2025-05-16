"""
Tests for percent-encoded characters in VCF files.
"""
import unittest
import tempfile
from pathlib import Path

from pygnome.parsers.vcf.vcf_reader import VcfReader
from pygnome.parsers.vcf.vcf_record import decode_percent_encoded


class TestPercentEncoding(unittest.TestCase):
    """Test cases for handling percent-encoded characters in VCF files."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary VCF file with percent-encoded characters
        self.temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.vcf')
        self.temp_file.write(b"""##fileformat=VCFv4.5
##INFO=<ID=DESC,Number=1,Type=String,Description="Description with special characters">
##INFO=<ID=TAGS,Number=.,Type=String,Description="Multiple tags">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t.\tA\tG\t30\tPASS\tDESC=Value%20with%20spaces%3B%20semicolons%3A%20colons%2C%20commas
1\t200\t.\tC\tT\t40\tPASS\tTAGS=tag1%2Ctag2,tag3%3Bwith%3Bsemicolons
1\t300\t.\tG\tA\t50\tPASS\tDESC=Value%25with%25percent%0Aand%0Dnewlines
""")
        self.temp_file.close()
        self.temp_path = Path(self.temp_file.name)
    
    def tearDown(self):
        """Clean up test fixtures."""
        # Remove the temporary file
        self.temp_path.unlink()
    
    def test_vcf_reader_with_percent_encoding(self):
        """Test that VcfReader correctly handles percent-encoded characters."""
        reader = VcfReader(self.temp_path)
        
        # Read all records
        records = list(reader)
        self.assertEqual(len(records), 3)
        
        # Check the first record
        self.assertEqual(records[0].get_pos(), 99)  # 0-based
        self.assertEqual(records[0].get_ref(), 'A')
        self.assertEqual(records[0].get_alt(), ['G'])
        self.assertEqual(records[0].get_info('DESC'), 'Value with spaces; semicolons: colons, commas')
        
        # Check the second record
        self.assertEqual(records[1].get_pos(), 199)  # 0-based
        self.assertEqual(records[1].get_ref(), 'C')
        self.assertEqual(records[1].get_alt(), ['T'])
        tags = records[1].get_info('TAGS')
        self.assertEqual(len(tags), 2)
        self.assertEqual(tags[0], 'tag1,tag2')  # This should be treated as a single value due to the escaped comma
        self.assertEqual(tags[1], 'tag3;with;semicolons')
        
        # Check the third record
        self.assertEqual(records[2].get_pos(), 299)  # 0-based
        self.assertEqual(records[2].get_ref(), 'G')
        self.assertEqual(records[2].get_alt(), ['A'])
        self.assertEqual(records[2].get_info('DESC'), 'Value%with%percent\nand\rnewlines')
    
    def test_decode_percent_encoded(self):
        """Test the decode_percent_encoded function."""
        # Test basic decoding
        self.assertEqual(decode_percent_encoded('Hello%20World'), 'Hello World')
        
        # Test decoding of special characters
        self.assertEqual(decode_percent_encoded('value%3Bwith%3Bsemicolons'), 'value;with;semicolons')
        self.assertEqual(decode_percent_encoded('value%3Awith%3Acolons'), 'value:with:colons')
        self.assertEqual(decode_percent_encoded('value%3Dwith%3Dequals'), 'value=with=equals')
        self.assertEqual(decode_percent_encoded('value%2Cwith%2Ccommas'), 'value,with,commas')
        self.assertEqual(decode_percent_encoded('value%25with%25percent'), 'value%with%percent')
        
        # Test decoding of newlines
        self.assertEqual(decode_percent_encoded('line1%0Aline2'), 'line1\nline2')
        self.assertEqual(decode_percent_encoded('line1%0Dline2'), 'line1\rline2')
        self.assertEqual(decode_percent_encoded('line1%0D%0Aline2'), 'line1\r\nline2')
        
        # Test decoding of tabs
        self.assertEqual(decode_percent_encoded('column1%09column2'), 'column1\tcolumn2')
        
        # Test with no encoded characters
        self.assertEqual(decode_percent_encoded('Hello World'), 'Hello World')
        
        # Test with invalid encoding (missing second hex digit)
        self.assertEqual(decode_percent_encoded('Invalid%2'), 'Invalid%2')
        
        # Test with invalid encoding (non-hex character)
        self.assertEqual(decode_percent_encoded('Invalid%ZZ'), 'Invalid%ZZ')


if __name__ == '__main__':
    unittest.main()