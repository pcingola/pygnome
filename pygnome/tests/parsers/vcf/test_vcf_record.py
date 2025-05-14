"""
Tests for the VcfRecord class.
"""
import unittest

from pygnome.parsers.vcf.vcf_header import VcfHeader
from pygnome.parsers.vcf.vcf_record import VcfRecord, Genotype


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
        self.assertEqual(self.record.get_chrom(), '20')
        self.assertEqual(self.record.get_pos(), 14369)  # 0-based
        self.assertEqual(self.record.get_vcf_pos(), 14370)  # 1-based
        self.assertEqual(self.record.get_id(), 'rs6054257')
        self.assertEqual(self.record.get_ref(), 'G')
        self.assertEqual(self.record.get_alt(), ['A'])
        self.assertEqual(self.record.get_qual(), 29.0)
        self.assertEqual(self.record.get_filter(), [])  # PASS means no filters
    
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
        self.assertEqual(self.record.get_format('GQ', 0), 48)
        self.assertEqual(self.record.get_format('DP', 0), 1)
        self.assertEqual(self.record.get_format('HQ', 0), [51, 51])
        
        # Check values for the second sample
        self.assertEqual(self.record.get_format('GQ', 1), 48)
        self.assertEqual(self.record.get_format('DP', 1), 8)
        self.assertEqual(self.record.get_format('HQ', 1), [51, 51])
        
        # Check values for the third sample
        self.assertEqual(self.record.get_format('GQ', 2), 43)
        self.assertEqual(self.record.get_format('DP', 2), 5)
        self.assertEqual(self.record.get_format('HQ', 2), [None, None])  # .,.
    
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
        self.assertEqual(self.multi_record.get_ref(), 'A')
        self.assertEqual(self.multi_record.get_alt(), ['G', 'T'])
        
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


if __name__ == '__main__':
    unittest.main()