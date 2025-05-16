"""
Tests for the VcfHeader class.
"""
import unittest
from pathlib import Path

from pygnome.parsers.vcf.vcf_header import VcfHeader


class TestVcfHeader(unittest.TestCase):
    """Test cases for the VcfHeader class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.header = VcfHeader()
        
        # Add some meta lines
        self.header.add_meta_line('##fileformat=VCFv4.5')
        self.header.add_meta_line('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">')
        self.header.add_meta_line('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">')
        self.header.add_meta_line('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">')
        self.header.add_meta_line('##FILTER=<ID=q10,Description="Quality below 10">')
        self.header.add_meta_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        self.header.add_meta_line('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">')
        self.header.add_meta_line('##contig=<ID=20,length=62435964,assembly=B36>')
        
        # Add sample names
        self.header.add_samples(['NA00001', 'NA00002', 'NA00003'])
    
    def test_file_format(self):
        """Test that the file format is correctly parsed."""
        self.assertEqual(self.header.file_format, 'VCFv4.5')
    
    def test_info_fields(self):
        """Test that INFO fields are correctly parsed."""
        # Check that we have the expected number of INFO fields
        self.assertEqual(len(self.header.info_fields), 3)
        
        # Check that the NS field is correctly parsed
        ns_field = self.header.get_info_field_definition('NS')
        if ns_field is None:
            raise ValueError("NS field is None")
        self.assertEqual(ns_field.id, 'NS')
        self.assertEqual(ns_field.number, '1')
        self.assertEqual(ns_field.type, 'Integer')
        self.assertEqual(ns_field.description, 'Number of Samples With Data')
        
        # Check that the AF field is correctly parsed
        af_field = self.header.get_info_field_definition('AF')
        if af_field is None:
            raise ValueError("AF field is None")
        self.assertEqual(af_field.id, 'AF')
        self.assertEqual(af_field.number, 'A')
        self.assertEqual(af_field.type, 'Float')
        self.assertEqual(af_field.description, 'Allele Frequency')
    
    def test_format_fields(self):
        """Test that FORMAT fields are correctly parsed."""
        # Check that we have the expected number of FORMAT fields
        self.assertEqual(len(self.header.format_fields), 2)
        
        # Check that the GT field is correctly parsed
        gt_field = self.header.get_format_field_definition('GT')
        if gt_field is None:
            raise ValueError("GT field is None")
        self.assertEqual(gt_field.id, 'GT')
        self.assertEqual(gt_field.number, '1')
        self.assertEqual(gt_field.type, 'String')
        self.assertEqual(gt_field.description, 'Genotype')
        
        # Check that the GQ field is correctly parsed
        gq_field = self.header.get_format_field_definition('GQ')
        if gq_field is None:
            raise ValueError("GQ field is None")
        self.assertEqual(gq_field.id, 'GQ')
        self.assertEqual(gq_field.number, '1')
        self.assertEqual(gq_field.type, 'Integer')
        self.assertEqual(gq_field.description, 'Genotype Quality')
    
    def test_filters(self):
        """Test that FILTER fields are correctly parsed."""
        # Check that we have the expected number of FILTER fields
        self.assertEqual(len(self.header.filters), 1)
        
        # Check that the q10 filter is correctly parsed
        q10_filter = self.header.filters.get('q10')
        if q10_filter is None:
            raise ValueError("q10 filter is None")
        self.assertEqual(q10_filter.id, 'q10')
        self.assertEqual(q10_filter.description, 'Quality below 10')
    
    def test_contigs(self):
        """Test that contig fields are correctly parsed."""
        # Check that we have the expected number of contigs
        self.assertEqual(len(self.header.contigs), 1)
        
        # Check that the contig is correctly parsed
        contig = self.header.contigs.get('20')
        if contig is None:
            raise ValueError("Contig 20 is None")
        self.assertEqual(contig.id, '20')
        self.assertEqual(contig.length, 62435964)
        self.assertEqual(contig.md5, None)
        
        # Check that the contigs list is correct
        self.assertEqual(self.header.get_contigs(), ['20'])
    
    def test_samples(self):
        """Test that sample names are correctly stored."""
        self.assertEqual(self.header.samples, ['NA00001', 'NA00002', 'NA00003'])
    
    def test_str(self):
        """Test the string representation of the header."""
        header_str = str(self.header)
        self.assertIn('##fileformat=VCFv4.5', header_str)
        self.assertIn('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">', header_str)
        self.assertIn('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003', header_str)


if __name__ == '__main__':
    unittest.main()