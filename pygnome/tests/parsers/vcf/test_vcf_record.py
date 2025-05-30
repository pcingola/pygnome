"""
Tests for the VcfRecord class.
"""
import unittest

from pygnome.parsers.vcf.vcf_header import VcfHeader
from pygnome.parsers.vcf.vcf_record import VcfRecord, Genotype
from pygnome.parsers.vcf.vcf_field_parser import decode_percent_encoded
from pygnome.genomics.variant import SNP, Insertion, Deletion, ComplexVariant


class TestVcfRecord(unittest.TestCase):
    """Test cases for the VcfRecord class."""
    
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
        self.header.add_meta_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        self.header.add_meta_line('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">')
        self.header.add_meta_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">')
        self.header.add_meta_line('##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">')
        
        # Add sample names
        self.header.add_samples(['NA00001', 'NA00002', 'NA00003'])
        
        # Create a record
        self.record_line = '20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.'
        self.record = VcfRecord(self.record_line, self.header)
        
        # Create a multi-allelic record
        self.multi_record_line = '20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4'
        self.multi_record = VcfRecord(self.multi_record_line, self.header)
    
    def test_fixed_fields(self):
        """Test that fixed fields are correctly parsed."""
        self.assertEqual(self.record._parse_chrom(), '20')
        self.assertEqual(self.record._parse_start(), 14369)  # 0-based
        self.assertEqual(self.record.pos, 14370)  # 1-based
        self.assertEqual(self.record._parse_id(), 'rs6054257')
        self.assertEqual(self.record._parse_ref(), 'G')
        self.assertEqual(self.record._parse_alt(), ['A'])
        self.assertEqual(self.record._parse_qual(), 29.0)
        self.assertEqual(self.record._parse_filter(), [])  # PASS means no filters
    
    def test_info_fields(self):
        """Test that INFO fields are correctly parsed."""
        self.assertTrue(self.record.has_info('NS'))
        self.assertEqual(self.record.get_info('NS'), 3)
        
        self.assertTrue(self.record.has_info('DP'))
        self.assertEqual(self.record.get_info('DP'), 14)
        
        self.assertTrue(self.record.has_info('AF'))
        self.assertEqual(self.record.get_info('AF'), 0.5)
        
        self.assertTrue(self.record.has_info('DB'))
        self.assertEqual(self.record.get_info('DB'), True)
        
        self.assertFalse(self.record.has_info('XX'))
        self.assertIsNone(self.record.get_info('XX'))
    
    def test_format_fields(self):
        """Test that FORMAT fields are correctly parsed."""
        self.assertTrue(self.record.has_format('GT'))
        self.assertTrue(self.record.has_format('GQ'))
        self.assertTrue(self.record.has_format('DP'))
        self.assertTrue(self.record.has_format('HQ'))
        
        # Check values for the first sample
        self.assertEqual(self.record.get_genotype_value('GQ', 0), 48)
        self.assertEqual(self.record.get_genotype_value('DP', 0), 1)
        self.assertEqual(self.record.get_genotype_value('HQ', 0), [51, 51])
        
        # Check values for the second sample
        self.assertEqual(self.record.get_genotype_value('GQ', 1), 48)
        self.assertEqual(self.record.get_genotype_value('DP', 1), 8)
        self.assertEqual(self.record.get_genotype_value('HQ', 1), [51, 51])
        
        # Check values for the third sample
        self.assertEqual(self.record.get_genotype_value('GQ', 2), 43)
        self.assertEqual(self.record.get_genotype_value('DP', 2), 5)
        self.assertEqual(self.record.get_genotype_value('HQ', 2), [None, None])  # .,.
    
    def test_genotypes(self):
        """Test that genotypes are correctly parsed."""
        genotypes = self.record.get_genotypes()
        
        # Check that we have the expected number of genotypes
        self.assertEqual(len(genotypes), 3)
        
        # Check the first genotype (0|0)
        self.assertEqual(genotypes[0].allele_indices, [0, 0])
        self.assertTrue(genotypes[0].phased)
        
        # Check the second genotype (1|0)
        self.assertEqual(genotypes[1].allele_indices, [1, 0])
        self.assertTrue(genotypes[1].phased)
        
        # Check the third genotype (1/1)
        self.assertEqual(genotypes[2].allele_indices, [1, 1])
        self.assertFalse(genotypes[2].phased)
    
    def test_multi_allelic(self):
        """Test parsing of multi-allelic records."""
        self.assertEqual(self.multi_record._parse_ref(), 'A')
        self.assertEqual(self.multi_record._parse_alt(), ['G', 'T'])
        
        # Check INFO fields with multiple values
        self.assertTrue(self.multi_record.has_info('AF'))
        self.assertEqual(self.multi_record.get_info('AF'), [0.333, 0.667])
        
        # Check genotypes
        genotypes = self.multi_record.get_genotypes()
        
        # Check the first genotype (1|2)
        self.assertEqual(genotypes[0].allele_indices, [1, 2])
        self.assertTrue(genotypes[0].phased)
        
        # Check the second genotype (2|1)
        self.assertEqual(genotypes[1].allele_indices, [2, 1])
        self.assertTrue(genotypes[1].phased)
        
        # Check the third genotype (2/2)
        self.assertEqual(genotypes[2].allele_indices, [2, 2])
        self.assertFalse(genotypes[2].phased)
    
    def test_alleles(self):
        """Test getting all alleles."""
        self.assertEqual(self.record.get_alleles(), ['G', 'A'])
        self.assertEqual(self.multi_record.get_alleles(), ['A', 'G', 'T'])
    
    def test_sample_names(self):
        """Test getting sample names."""
        self.assertEqual(self.record.get_sample_names(), ['NA00001', 'NA00002', 'NA00003'])
    
    def test_str(self):
        """Test the string representation of the record."""
        self.assertEqual(str(self.record), self.record_line)
    
    def test_iterator(self):
        """Test the iterator functionality of VcfRecord."""
        # Test with SNP record
        variants = list(self.record)
        self.assertEqual(len(variants), 1)
        self.assertIsInstance(variants[0], SNP)
        self.assertEqual(variants[0].ref, "G")
        self.assertEqual(variants[0].alt, "A")  # Should be a string, not a list
        
        # Test with multi-allelic record
        variants = list(self.multi_record)
        self.assertEqual(len(variants), 2)
        
        # First variant should be an insertion (A -> G)
        self.assertEqual(variants[0].ref, "A")
        self.assertEqual(variants[0].alt, "G")  # Should be a string, not a list
        
        # Second variant should be an insertion (A -> T)
        self.assertEqual(variants[1].ref, "A")
        self.assertEqual(variants[1].alt, "T")  # Should be a string, not a list
    
    def test_percent_encoded_characters(self):
        """Test handling of percent-encoded characters in INFO fields."""
        # Test the decode_percent_encoded function
        self.assertEqual(decode_percent_encoded("Hello%20World"), "Hello World")
        self.assertEqual(decode_percent_encoded("value%3Bwith%3Bsemicolons"), "value;with;semicolons")
        self.assertEqual(decode_percent_encoded("value%3Awith%3Acolons"), "value:with:colons")
        self.assertEqual(decode_percent_encoded("value%3Dwith%3Dequals"), "value=with=equals")
        self.assertEqual(decode_percent_encoded("value%2Cwith%2Ccommas"), "value,with,commas")
        self.assertEqual(decode_percent_encoded("value%25with%25percent"), "value%with%percent")
        
        # Create a header with an INFO field that might contain escaped characters
        header = VcfHeader()
        header.add_meta_line('##fileformat=VCFv4.5')
        header.add_meta_line('##INFO=<ID=DESC,Number=1,Type=String,Description="Description with special characters">')
        
        # Create a record with percent-encoded characters in the INFO field
        record_line = '1\t100\t.\tA\tG\t30\tPASS\tDESC=Value%20with%20spaces%3B%20semicolons%3A%20colons%2C%20commas\tGT\t0/1'
        record = VcfRecord(record_line, header)
        
        # Check that the INFO field is correctly decoded
        self.assertEqual(record.get_info('DESC'), "Value with spaces; semicolons: colons, commas")
        
        # Test with multiple values separated by commas
        header.add_meta_line('##INFO=<ID=TAGS,Number=.,Type=String,Description="Multiple tags">')
        record_line = '1\t100\t.\tA\tG\t30\tPASS\tTAGS=tag1%2Ctag2,tag3%3Bwith%3Bsemicolons\tGT\t0/1'
        record = VcfRecord(record_line, header)
        
        # Check that the comma-separated values are correctly parsed
        tags = record.get_info('TAGS')
        self.assertEqual(len(tags), 2)
        self.assertEqual(tags[0], "tag1,tag2")  # This should be treated as a single value due to the escaped comma
        self.assertEqual(tags[1], "tag3;with;semicolons")


    def test_lazy_genotype_parsing(self):
        """Test that genotype parsing is lazy and only happens when needed."""
        # Create a VCF record with many samples
        header = VcfHeader()
        header.add_meta_line('##fileformat=VCFv4.5')
        header.add_meta_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        header.add_meta_line('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">')
        header.add_meta_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">')
        
        # Add sample names
        sample_names = []
        for i in range(100):
            sample_names.append(f"SAMPLE{i}")
        header.add_samples(sample_names)
        
        # Create a VCF line with many samples
        line = "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP"
        for i in range(100):
            line += f"\t0/1:{i}:10"
        
        # Create a VCF record
        record = VcfRecord(line, header)
        
        # Verify that the record has genotypes
        self.assertTrue(record.genotypes.has_genotypes)
        
        # Verify that the fields haven't been parsed yet in the genotypes object
        self.assertIsNone(record.genotypes._fields)
        
        # Access a specific genotype value - this should trigger parsing only for that sample
        gq_value = record.get_genotype_value("GQ", 50)
        
        # Verify the value
        self.assertEqual(gq_value, 50)
        
        # Verify that only one sample was parsed
        self.assertEqual(len(record.genotypes._sample_data_cache), 1)
        
        # Access another sample
        gq_value = record.get_genotype_value("GQ", 75)
        self.assertEqual(gq_value, 75)
        
        # Verify that only two samples were parsed
        self.assertEqual(len(record.genotypes._sample_data_cache), 2)
        
        # Get all genotypes - this should parse all samples
        all_genotypes = record.get_genotypes()
        
        # Verify that we have the expected number of genotypes
        self.assertEqual(len(all_genotypes), 100)
        
        # Verify that all samples have been parsed
        self.assertEqual(len(record.genotypes._sample_data_cache), 100)
        
    def test_lazy_string_conversion(self):
        """Test that string conversion is lazy and uses the cached raw line if no changes."""
        # Create a VCF record
        line = "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t0/1:50:10\t1/1:60:20"
        record = VcfRecord(line, self.header)
        
        # Convert to string without making any changes
        record_str = str(record)
        
        # Verify that the original line is returned
        self.assertEqual(record_str, line)
        
        # Make a change to the genotype
        record.set_genotype_value("GQ", 55, 0)
        
        # Convert to string again
        record_str = str(record)
        
        # Verify that the string has been updated
        self.assertNotEqual(record_str, line)
        self.assertIn("0/1:55:10", record_str)


if __name__ == '__main__':
    unittest.main()