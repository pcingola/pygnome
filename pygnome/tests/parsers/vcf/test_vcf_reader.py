"""
Tests for the VcfReader class.
"""
import os
import unittest
from pathlib import Path

from pygnome.parsers.vcf.vcf_reader import VcfReader
from pygnome.parsers.vcf.vcf_record import VcfRecord


class TestVcfReader(unittest.TestCase):
    """Test cases for the VcfReader class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Get the path to the sample VCF file
        self.sample_vcf_path = Path(__file__).parent / "data" / "sample.vcf"
        
        # Make sure the file exists
        self.assertTrue(self.sample_vcf_path.exists(), f"Sample VCF file not found: {self.sample_vcf_path}")
    
    def test_header_parsing(self):
        """Test that the header is correctly parsed."""
        with VcfReader(self.sample_vcf_path) as reader:
            # Check file format
            self.assertEqual(reader.header.file_format, "VCFv4.5")
            
            # Check INFO fields
            self.assertEqual(len(reader.header.info_fields), 6)
            self.assertIn("NS", reader.header.info_fields)
            self.assertIn("DP", reader.header.info_fields)
            self.assertIn("AF", reader.header.info_fields)
            self.assertIn("AA", reader.header.info_fields)
            self.assertIn("DB", reader.header.info_fields)
            self.assertIn("H2", reader.header.info_fields)
            
            # Check FORMAT fields
            self.assertEqual(len(reader.header.format_fields), 4)
            self.assertIn("GT", reader.header.format_fields)
            self.assertIn("GQ", reader.header.format_fields)
            self.assertIn("DP", reader.header.format_fields)
            self.assertIn("HQ", reader.header.format_fields)
            
            # Check FILTER fields
            self.assertEqual(len(reader.header.filters), 2)
            self.assertIn("q10", reader.header.filters)
            self.assertIn("s50", reader.header.filters)
            
            # Check contigs
            self.assertEqual(len(reader.header.contigs), 1)
            self.assertIn("20", reader.header.contigs)
            
            # Check samples
            self.assertEqual(reader.header.samples, ["NA00001", "NA00002", "NA00003"])
    
    def test_record_iteration(self):
        """Test iterating through records."""
        with VcfReader(self.sample_vcf_path) as reader:
            records = list(reader)
            
            # Check that we have the expected number of records
            self.assertEqual(len(records), 5)
            
            # Check that all records are VcfRecord objects
            for record in records:
                self.assertIsInstance(record, VcfRecord)
            
            # Check the first record
            first_record = records[0]
            self.assertEqual(first_record._parse_chrom(), "20")
            self.assertEqual(first_record._parse_start(), 14369)  # 0-based
            self.assertEqual(first_record._parse_id(), "rs6054257")
            self.assertEqual(first_record._parse_ref(), "G")
            self.assertEqual(first_record._parse_alt(), ["A"])
            
            # Check the last record
            last_record = records[-1]
            self.assertEqual(last_record._parse_chrom(), "20")
            self.assertEqual(last_record._parse_start(), 1234566)  # 0-based
            self.assertEqual(last_record._parse_id(), "microsat1")
            self.assertEqual(last_record._parse_ref(), "GTC")
            self.assertEqual(last_record._parse_alt(), ["G", "GTCT"])
    
    def test_fetch(self):
        """Test fetching records by region."""
        with VcfReader(self.sample_vcf_path) as reader:
            # Fetch records in a region
            records = list(reader.fetch("20", 14000, 15000))
            
            # Check that we have the expected number of records
            self.assertEqual(len(records), 1)
            
            # Check the record
            record = records[0]
            self.assertEqual(record._parse_chrom(), "20")
            self.assertEqual(record._parse_start(), 14369)  # 0-based
            self.assertEqual(record._parse_id(), "rs6054257")
            
            # Fetch records in another region
            records = list(reader.fetch("20", 1000000, 2000000))
            
            # Check that we have the expected number of records
            self.assertEqual(len(records), 3)
            
            # Check the records
            self.assertEqual(records[0]._parse_start(), 1110695)  # 0-based
            self.assertEqual(records[1]._parse_start(), 1230236)  # 0-based
            self.assertEqual(records[2]._parse_start(), 1234566)  # 0-based
    

    def test_context_manager(self):
        """Test using the reader as a context manager."""
        # Open the reader using a context manager
        with VcfReader(self.sample_vcf_path) as reader:
            # Check that the reader is initialized
            self.assertIsNotNone(reader.header)
            self.assertEqual(reader.header.file_format, "VCFv4.5")
            
            # Read some records
            records = list(reader)
            self.assertEqual(len(records), 5)
        
        # The file should be closed after the context manager exits
        # We can't directly test this, but we can check that we can open it again
        with VcfReader(self.sample_vcf_path) as reader:
            records = list(reader)
            self.assertEqual(len(records), 5)
    
    def test_get_samples(self):
        """Test getting sample names."""
        with VcfReader(self.sample_vcf_path) as reader:
            self.assertEqual(reader.get_samples(), ["NA00001", "NA00002", "NA00003"])
    
    def test_get_contigs(self):
        """Test getting contig names."""
        with VcfReader(self.sample_vcf_path) as reader:
            self.assertEqual(reader.get_contigs(), ["20"])


if __name__ == "__main__":
    unittest.main()