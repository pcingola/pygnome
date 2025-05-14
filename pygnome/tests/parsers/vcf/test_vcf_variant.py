"""
Tests for the VcfVariant class.
"""
import unittest

from pygnome.parsers.vcf.vcf_header import VcfHeader
from pygnome.parsers.vcf.vcf_record import VcfRecord
from pygnome.parsers.vcf.vcf_variant import VcfVariant, VariantType


class TestVcfVariant(unittest.TestCase):
    """Test cases for the VcfVariant class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a header
        self.header = VcfHeader()
        
        # Add meta lines
        self.header.add_meta_line('##fileformat=VCFv4.5')
        self.header.add_meta_line('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">')
        self.header.add_meta_line('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">')
        self.header.add_meta_line('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">')
        self.header.add_meta_line('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">')
        self.header.add_meta_line('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">')
        self.header.add_meta_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        
        # Add sample names
        self.header.add_samples(['NA00001', 'NA00002'])
        
        # Create records for different variant types
        
        # SNP
        self.snp_line = '20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5\tGT\t0|0\t1|0'
        self.snp_record = VcfRecord(self.snp_line, self.header)
        self.snp_variant = VcfVariant(self.snp_record)
        
        # Insertion
        self.ins_line = '20\t17330\t.\tT\tTA\t3\tPASS\tNS=3;DP=11;AF=0.017\tGT\t0|0\t0|1'
        self.ins_record = VcfRecord(self.ins_line, self.header)
        self.ins_variant = VcfVariant(self.ins_record)
        
        # Deletion
        self.del_line = '20\t1110696\trs6040355\tAT\tA\t67\tPASS\tNS=2;DP=10;AF=0.333\tGT\t1|0\t0|1'
        self.del_record = VcfRecord(self.del_line, self.header)
        self.del_variant = VcfVariant(self.del_record)
        
        # Multi-allelic
        self.multi_line = '20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9\tGT\t0/1\t0/2'
        self.multi_record = VcfRecord(self.multi_line, self.header)
        self.multi_variant = VcfVariant(self.multi_record)
        
        # Structural variant
        self.sv_line = '20\t1234568\tsv1\tG\t<DEL>\t50\tPASS\tSVLEN=1000\tGT\t0/1\t0/1'
        self.sv_record = VcfRecord(self.sv_line, self.header)
        self.sv_variant = VcfVariant(self.sv_record)
        
        # Breakend
        self.bnd_line = '20\t1234569\tbnd1\tG\tG[2:321682[\t50\tPASS\tNS=3;DP=9\tGT\t0/1\t0/1'
        self.bnd_record = VcfRecord(self.bnd_line, self.header)
        self.bnd_variant = VcfVariant(self.bnd_record)
    
    def test_basic_properties(self):
        """Test basic properties of variants."""
        # SNP
        self.assertEqual(self.snp_variant.get_chrom(), '20')
        self.assertEqual(self.snp_variant.get_pos(), 14369)  # 0-based
        self.assertEqual(self.snp_variant.get_end(), 14370)  # 0-based, exclusive
        self.assertEqual(self.snp_variant.get_id(), 'rs6054257')
        self.assertEqual(self.snp_variant.get_ref(), 'G')
        self.assertEqual(self.snp_variant.get_alt(), ['A'])
        self.assertEqual(self.snp_variant.get_alleles(), ['G', 'A'])
        
        # Insertion
        self.assertEqual(self.ins_variant.get_chrom(), '20')
        self.assertEqual(self.ins_variant.get_pos(), 17329)  # 0-based
        self.assertEqual(self.ins_variant.get_end(), 17330)  # 0-based, exclusive
        self.assertEqual(self.ins_variant.get_ref(), 'T')
        self.assertEqual(self.ins_variant.get_alt(), ['TA'])
        
        # Deletion
        self.assertEqual(self.del_variant.get_chrom(), '20')
        self.assertEqual(self.del_variant.get_pos(), 1110695)  # 0-based
        self.assertEqual(self.del_variant.get_end(), 1110697)  # 0-based, exclusive
        self.assertEqual(self.del_variant.get_ref(), 'AT')
        self.assertEqual(self.del_variant.get_alt(), ['A'])
    
    def test_variant_types(self):
        """Test variant type detection."""
        # SNP
        self.assertTrue(self.snp_variant.is_snp())
        self.assertFalse(self.snp_variant.is_indel())
        self.assertFalse(self.snp_variant.is_insertion())
        self.assertFalse(self.snp_variant.is_deletion())
        self.assertFalse(self.snp_variant.is_structural_variant())
        self.assertFalse(self.snp_variant.is_breakend())
        self.assertFalse(self.snp_variant.is_multi_allelic())
        self.assertEqual(self.snp_variant.get_variant_type(), VariantType.SNP)
        
        # Insertion
        self.assertFalse(self.ins_variant.is_snp())
        self.assertTrue(self.ins_variant.is_indel())
        self.assertTrue(self.ins_variant.is_insertion())
        self.assertFalse(self.ins_variant.is_deletion())
        self.assertFalse(self.ins_variant.is_structural_variant())
        self.assertFalse(self.ins_variant.is_breakend())
        self.assertFalse(self.ins_variant.is_multi_allelic())
        self.assertEqual(self.ins_variant.get_variant_type(), VariantType.INS)
        
        # Deletion
        self.assertFalse(self.del_variant.is_snp())
        self.assertTrue(self.del_variant.is_indel())
        self.assertFalse(self.del_variant.is_insertion())
        self.assertTrue(self.del_variant.is_deletion())
        self.assertFalse(self.del_variant.is_structural_variant())
        self.assertFalse(self.del_variant.is_breakend())
        self.assertFalse(self.del_variant.is_multi_allelic())
        self.assertEqual(self.del_variant.get_variant_type(), VariantType.DEL)
        
        # Multi-allelic
        self.assertFalse(self.multi_variant.is_snp())
        self.assertTrue(self.multi_variant.is_indel())
        self.assertTrue(self.multi_variant.is_insertion())
        self.assertTrue(self.multi_variant.is_deletion())
        self.assertFalse(self.multi_variant.is_structural_variant())
        self.assertFalse(self.multi_variant.is_breakend())
        self.assertTrue(self.multi_variant.is_multi_allelic())
        
        # Structural variant
        self.assertFalse(self.sv_variant.is_snp())
        self.assertFalse(self.sv_variant.is_indel())
        self.assertFalse(self.sv_variant.is_insertion())
        self.assertFalse(self.sv_variant.is_deletion())
        self.assertTrue(self.sv_variant.is_structural_variant())
        self.assertFalse(self.sv_variant.is_breakend())
        self.assertFalse(self.sv_variant.is_multi_allelic())
        self.assertEqual(self.sv_variant.get_variant_type(), "DEL")  # From the symbolic allele
        
        # Breakend
        self.assertFalse(self.bnd_variant.is_snp())
        self.assertFalse(self.bnd_variant.is_indel())
        self.assertFalse(self.bnd_variant.is_insertion())
        self.assertFalse(self.bnd_variant.is_deletion())
        self.assertFalse(self.bnd_variant.is_structural_variant())
        self.assertTrue(self.bnd_variant.is_breakend())
        self.assertFalse(self.bnd_variant.is_multi_allelic())
        self.assertEqual(self.bnd_variant.get_variant_type(), VariantType.BND)
    
    def test_genotypes(self):
        """Test getting genotypes."""
        # SNP
        genotypes = self.snp_variant.get_genotypes()
        self.assertEqual(len(genotypes), 2)
        self.assertEqual(genotypes[0].allele_indices, [0, 0])
        self.assertTrue(genotypes[0].phased)
        self.assertEqual(genotypes[1].allele_indices, [1, 0])
        self.assertTrue(genotypes[1].phased)
    
    def test_sample_names(self):
        """Test getting sample names."""
        self.assertEqual(self.snp_variant.get_sample_names(), ['NA00001', 'NA00002'])
    
    def test_str(self):
        """Test the string representation of the variant."""
        self.assertEqual(str(self.snp_variant), '20:14370:G>A')
        self.assertEqual(str(self.ins_variant), '20:17330:T>TA')
        self.assertEqual(str(self.del_variant), '20:1110696:AT>A')
        self.assertEqual(str(self.multi_variant), '20:1234567:GTC>G,GTCT')


if __name__ == '__main__':
    unittest.main()