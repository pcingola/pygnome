"""
Tests for adding and modifying INFO fields in VCF records.
"""
import unittest

from pygnome.parsers.vcf.vcf_header import VcfHeader, FieldType
from pygnome.parsers.vcf.vcf_record import VcfRecord
from pygnome.parsers.vcf.vcf_field_parser import encode_percent_encoded


class TestInfoModification(unittest.TestCase):
    """Test cases for adding and modifying INFO fields in VCF records."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a header
        self.header = VcfHeader()
        
        # Add meta lines
        self.header.add_meta_line('##fileformat=VCFv4.5')
        self.header.add_meta_line('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">')
        self.header.add_meta_line('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">')
        self.header.add_meta_line('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">')
        self.header.add_meta_line('##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">')
        self.header.add_meta_line('##INFO=<ID=DESC,Number=1,Type=String,Description="Description field">')
        self.header.add_meta_line('##INFO=<ID=TAGS,Number=.,Type=String,Description="Multiple tags">')
        self.header.add_meta_line('##INFO=<ID=TEST,Number=1,Type=Integer,Description="Test field for unit tests">')
        self.header.add_meta_line('##INFO=<ID=TEST_FLAG,Number=0,Type=Flag,Description="Test flag field for unit tests">')
        self.header.add_meta_line('##INFO=<ID=TEST_FLAG_FALSE,Number=0,Type=Flag,Description="Test flag field for unit tests">')
        self.header.add_meta_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        
        # Add sample names
        self.header.add_samples(['NA00001'])
        
        # Create a record with some INFO fields
        self.record_line = '20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB\tGT\t0/1'
        self.record = VcfRecord(self.record_line, self.header)
        
        # Create a record with no INFO fields
        self.empty_record_line = '20\t14371\t.\tG\tA\t29\tPASS\t.\tGT\t0/1'
        self.empty_record = VcfRecord(self.empty_record_line, self.header)
    
    def test_encode_percent_encoded(self):
        """Test the encode_percent_encoded function."""
        # Test basic encoding
        self.assertEqual(encode_percent_encoded('Hello World'), 'Hello%20World')
        
        # Test encoding of special characters
        self.assertEqual(encode_percent_encoded('value;with;semicolons'), 'value%3Bwith%3Bsemicolons')
        self.assertEqual(encode_percent_encoded('value:with:colons'), 'value%3Awith%3Acolons')
        self.assertEqual(encode_percent_encoded('value=with=equals'), 'value%3Dwith%3Dequals')
        self.assertEqual(encode_percent_encoded('value,with,commas'), 'value%2Cwith%2Ccommas')
        self.assertEqual(encode_percent_encoded('value%with%percent'), 'value%25with%25percent')
        
        # Test encoding of newlines
        self.assertEqual(encode_percent_encoded('line1\nline2'), 'line1%0Aline2')
        self.assertEqual(encode_percent_encoded('line1\rline2'), 'line1%0Dline2')
        self.assertEqual(encode_percent_encoded('line1\r\nline2'), 'line1%0D%0Aline2')
        
        # Test encoding of tabs
        self.assertEqual(encode_percent_encoded('column1\tcolumn2'), 'column1%09column2')
        
        # Test with None - this should be handled in _format_info_value, not in encode_percent_encoded
        # Removing this test case as encode_percent_encoded expects a string
    
    def test_add_info_to_existing_record(self):
        """Test adding a new INFO field to a record with existing INFO fields."""
        # Add a new integer field
        self.record.add_info('TEST', 42)
        
        # Check that the field was added
        self.assertTrue(self.record.has_info('TEST'))
        self.assertEqual(self.record.get_info('TEST'), 42)
        
        # Convert to string to update the raw line
        record_str = str(self.record)
        
        # Check that the raw line was updated
        self.assertIn('TEST=42', record_str)
        
        # Check that existing fields are preserved
        self.assertEqual(self.record.get_info('NS'), 3)
        self.assertEqual(self.record.get_info('DP'), 14)
        self.assertEqual(self.record.get_info('AF'), 0.5)
        self.assertTrue(self.record.get_info('DB'))
    
    def test_add_info_to_empty_record(self):
        """Test adding a new INFO field to a record with no INFO fields."""
        # Add a new integer field
        self.empty_record.add_info('NS', 3)
        
        # Check that the field was added
        self.assertTrue(self.empty_record.has_info('NS'))
        self.assertEqual(self.empty_record.get_info('NS'), 3)
        
        # Convert to string to update the fields
        str(self.empty_record)
        
        # Check that the fields were updated
        self.assertEqual(self.empty_record._fields[7], 'NS=3')
        
        # Convert to string to update the raw line
        record_str = str(self.empty_record)
        
        # Check that the raw line was updated
        self.assertIn('NS=3', record_str)
    
    def test_add_flag_info(self):
        """Test adding a flag INFO field."""
        # Add a new flag field
        self.record.add_info('TEST_FLAG', True)
        
        # Check that the field was added
        self.assertTrue(self.record.has_info('TEST_FLAG'))
        self.assertTrue(self.record.get_info('TEST_FLAG'))
        
        # Convert to string to update the raw line
        record_str = str(self.record)
        
        # Check that the raw line was updated
        self.assertIn('TEST_FLAG', record_str)
        
        # Add a flag field with False value (should not be added)
        self.record.add_info('TEST_FLAG_FALSE', False)
        
        # Check that the field was not added
        self.assertFalse(self.record.has_info('TEST_FLAG_FALSE'))
    
    def test_modify_existing_info(self):
        """Test modifying an existing INFO field."""
        # Modify an existing integer field
        self.record.set_info('NS', 5)
        
        # Check that the field was updated
        self.assertEqual(self.record.get_info('NS'), 5)
        
        # Convert to string to update the raw line
        record_str = str(self.record)
        
        # Check that the raw line was updated
        self.assertIn('NS=5', record_str)
        
        # Check that other fields are preserved
        self.assertEqual(self.record.get_info('DP'), 14)
        self.assertEqual(self.record.get_info('AF'), 0.5)
        self.assertTrue(self.record.get_info('DB'))
    
    def test_add_string_with_special_chars(self):
        """Test adding a string with special characters that need encoding."""
        # Add a string with special characters
        special_string = "Value with spaces; semicolons: colons, commas"
        self.record.add_info('DESC', special_string)
        
        # Check that the field was added
        self.assertTrue(self.record.has_info('DESC'))
        self.assertEqual(self.record.get_info('DESC'), special_string)
        
        # Check that the raw line contains the encoded string
        encoded_string = encode_percent_encoded(special_string)
        
        # Convert to string to update the raw line
        record_str = str(self.record)
        
        # Check that the raw line contains the encoded string
        self.assertIn(f'DESC={encoded_string}', record_str)
    
    def test_add_list_of_values(self):
        """Test adding a list of values."""
        # Add a list of floats
        self.record.add_info('AF', [0.1, 0.2, 0.3])
        
        # Check that the field was updated
        self.assertEqual(self.record.get_info('AF'), [0.1, 0.2, 0.3])
        
        # Convert to string to update the raw line
        record_str = str(self.record)
        
        # Check that the raw line was updated
        self.assertIn('AF=0.1,0.2,0.3', record_str)
    
    def test_add_list_with_special_chars(self):
        """Test adding a list of strings with special characters."""
        # Add a list of strings with special characters
        tags = ["tag1,tag2", "tag3;with;semicolons"]
        self.record.add_info('TAGS', tags)
        
        # Check that the field was added
        self.assertTrue(self.record.has_info('TAGS'))
        self.assertEqual(self.record.get_info('TAGS'), tags)
        
        # Check that the raw line contains the encoded strings
        encoded_tags = f"TAGS={encode_percent_encoded(tags[0])},{encode_percent_encoded(tags[1])}"
        record_str = str(self.record)
        self.assertIn(encoded_tags, record_str)


if __name__ == '__main__':
    unittest.main()