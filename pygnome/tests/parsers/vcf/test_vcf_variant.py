"""
Tests for variants created from VcfRecord.
"""
import unittest

from pygnome.parsers.vcf.vcf_header import VcfHeader
from pygnome.parsers.vcf.vcf_record import VcfRecord
from pygnome.genomics.variant import (
    VariantType, SNP, Insertion, Deletion,
    Duplication, Inversion, Translocation, ComplexVariant
)


class TestVcfVariant(unittest.TestCase):
    """Test cases for variants created from VcfRecord."""
    
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
        self.header.add_meta_line('##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate">')
        self.header.add_meta_line('##INFO=<ID=POS2,Number=1,Type=Integer,Description="Position for END coordinate">')
        self.header.add_meta_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        
        # Add sample names
        self.header.add_samples(['NA00001', 'NA00002'])
        
        # Define VCF record lines for different variant types
        
        # SNP
        self.snp_line = '20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5\tGT\t0|0\t1|0'
        self.snp_record = VcfRecord(self.snp_line, self.header)
        
        # Insertion
        self.ins_line = '20\t17330\t.\tT\tTA\t3\tPASS\tNS=3;DP=11;AF=0.017\tGT\t0|0\t0|1'
        self.ins_record = VcfRecord(self.ins_line, self.header)
        
        # Deletion
        self.del_line = '20\t1110696\trs6040355\tAT\tA\t67\tPASS\tNS=2;DP=10;AF=0.333\tGT\t1|0\t0|1'
        self.del_record = VcfRecord(self.del_line, self.header)
        
        # Multi-allelic
        self.multi_line = '20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9\tGT\t0/1\t0/2'
        self.multi_record = VcfRecord(self.multi_line, self.header)
        
        # Structural variants
        self.del_sv_line = '20\t1234568\tsv1\tG\t<DEL>\t50\tPASS\tSVLEN=1000\tGT\t0/1\t0/1'
        self.del_sv_record = VcfRecord(self.del_sv_line, self.header)
        
        self.dup_sv_line = '20\t1234570\tsv2\tG\t<DUP>\t50\tPASS\tSVLEN=500\tGT\t0/1\t0/1'
        self.dup_sv_record = VcfRecord(self.dup_sv_line, self.header)
        
        self.inv_sv_line = '20\t1234571\tsv3\tG\t<INV>\t50\tPASS\tSVLEN=300\tGT\t0/1\t0/1'
        self.inv_sv_record = VcfRecord(self.inv_sv_line, self.header)
        
        self.tra_sv_line = '20\t1234572\tsv4\tG\t<TRA>\t50\tPASS\tCHR2=X;POS2=15000\tGT\t0/1\t0/1'
        self.tra_sv_record = VcfRecord(self.tra_sv_line, self.header)
        
        # Breakend
        self.bnd_line = '20\t1234569\tbnd1\tG\tG[2:321682[\t50\tPASS\tNS=3;DP=9\tGT\t0/1\t0/1'
        self.bnd_record = VcfRecord(self.bnd_line, self.header)
    
    def test_snp_variant_creation(self):
        """Test creation of SNP variant from VcfRecord."""
        # Get the SNP variant using the iterator
        variants = list(self.snp_record)
        self.assertEqual(len(variants), 1)
        
        snp_variant = variants[0]
        
        # Verify it's the correct type
        self.assertIsInstance(snp_variant, SNP)
        
        # Test basic properties
        self.assertEqual(snp_variant.chrom, '20')
        self.assertEqual(snp_variant.start, 14369)  # 0-based
        self.assertEqual(snp_variant.end, 14370)  # 0-based, exclusive
        self.assertEqual(snp_variant.id, 'rs6054257')
        self.assertEqual(snp_variant.ref, 'G')
        self.assertEqual(snp_variant.alt, 'A')


    def test_insertion_variant_creation(self):
        """Test creation of Insertion variant from VcfRecord."""
        # Get the Insertion variant using the iterator
        variants = list(self.ins_record)
        self.assertEqual(len(variants), 1)
        
        ins_variant = variants[0]
        
        # Verify it's the correct type
        self.assertIsInstance(ins_variant, Insertion)
        
        # Test basic properties
        self.assertEqual(ins_variant.chrom, '20')
        self.assertEqual(ins_variant.start, 17329)  # 0-based
        self.assertEqual(ins_variant.end, 17330)  # 0-based, exclusive
        self.assertEqual(ins_variant.ref, 'T')
        self.assertEqual(ins_variant.alt, 'TA')
        
        # Test insertion-specific property
        self.assertIsInstance(ins_variant, Insertion)  # Ensure correct type before accessing specific property
        if isinstance(ins_variant, Insertion):
            self.assertEqual(ins_variant.inserted_sequence, 'A')
    
    def test_deletion_variant_creation(self):
        """Test creation of Deletion variant from VcfRecord."""
        # Get the Deletion variant using the iterator
        variants = list(self.del_record)
        self.assertEqual(len(variants), 1)
        
        del_variant = variants[0]
        
        # Verify it's the correct type
        self.assertIsInstance(del_variant, Deletion)
        
        # Test basic properties
        self.assertEqual(del_variant.chrom, '20')
        self.assertEqual(del_variant.start, 1110695)  # 0-based
        self.assertEqual(del_variant.end, 1110697)  # 0-based, exclusive
        self.assertEqual(del_variant.ref, 'AT')
        self.assertEqual(del_variant.alt, 'A')
        
        # Test deletion-specific property
        self.assertIsInstance(del_variant, Deletion)  # Ensure correct type before accessing specific property
        if isinstance(del_variant, Deletion):
            self.assertEqual(del_variant.deleted_sequence, 'T')
    
    def test_multi_allelic_variant_creation(self):
        """Test creation of variants from multi-allelic VcfRecord."""
        # Get all variants using the iterator
        variants = list(self.multi_record)
        self.assertEqual(len(variants), 2)
        
        # First variant (deletion)
        del_variant = variants[0]
        self.assertIsInstance(del_variant, Deletion)
        self.assertEqual(del_variant.ref, 'GTC')
        self.assertEqual(del_variant.alt, 'G')
        
        # Second variant (insertion)
        ins_variant = variants[1]
        self.assertIsInstance(ins_variant, Insertion)
        self.assertEqual(ins_variant.ref, 'GTC')
        self.assertEqual(ins_variant.alt, 'GTCT')
        self.assertIsInstance(ins_variant, Insertion)  # Ensure correct type before accessing specific property
        if isinstance(ins_variant, Insertion):
            self.assertEqual(ins_variant.inserted_sequence, 'T')
    
    def test_structural_variant_creation(self):
        """Test creation of structural variants from VcfRecord."""
        # Test deletion structural variant
        del_variants = list(self.del_sv_record)
        self.assertEqual(len(del_variants), 1)
        del_variant = del_variants[0]
        
        # DEL structural variants should be represented as Deletion
        self.assertIsInstance(del_variant, Deletion)
        self.assertEqual(del_variant.chrom, '20')
        self.assertEqual(del_variant.ref, 'G')
        self.assertEqual(del_variant.alt, '<DEL>')
        
        # Test duplication structural variant
        dup_variants = list(self.dup_sv_record)
        self.assertEqual(len(dup_variants), 1)
        dup_variant = dup_variants[0]
        
        # DUP structural variants should be represented as Duplication
        self.assertIsInstance(dup_variant, Duplication)
        self.assertEqual(dup_variant.chrom, '20')
        self.assertEqual(dup_variant.ref, 'G')
        self.assertEqual(dup_variant.alt, '<DUP>')
        # Check dup_length only if it's a Duplication
        if isinstance(dup_variant, Duplication):
            self.assertEqual(dup_variant.dup_length, 500)
        
        # Test inversion structural variant
        inv_variants = list(self.inv_sv_record)
        self.assertEqual(len(inv_variants), 1)
        inv_variant = inv_variants[0]
        
        # INV structural variants should be represented as Inversion
        self.assertIsInstance(inv_variant, Inversion)
        self.assertEqual(inv_variant.chrom, '20')
        self.assertEqual(inv_variant.ref, 'G')
        self.assertEqual(inv_variant.alt, '<INV>')
        # Check inv_length only if it's an Inversion
        if isinstance(inv_variant, Inversion):
            self.assertEqual(inv_variant.inv_length, 300)
        
        # Test translocation structural variant
        tra_variants = list(self.tra_sv_record)
        self.assertEqual(len(tra_variants), 1)
        tra_variant = tra_variants[0]
        
        # TRA structural variants should be represented as Translocation
        self.assertIsInstance(tra_variant, Translocation)
        self.assertEqual(tra_variant.chrom, '20')
        self.assertEqual(tra_variant.ref, 'G')
        self.assertEqual(tra_variant.alt, '<TRA>')
        # Check destination properties only if it's a Translocation
        if isinstance(tra_variant, Translocation):
            self.assertEqual(tra_variant.dest_chrom, 'X')
            self.assertEqual(tra_variant.dest_pos, 15000)
    
    def test_breakend_variant_creation(self):
        """Test creation of breakend variant from VcfRecord."""
        # Get the breakend variant using the iterator
        variants = list(self.bnd_record)
        self.assertEqual(len(variants), 1)
        
        bnd_variant = variants[0]
        
        # Breakend variants should be represented as Translocation when possible
        self.assertIsInstance(bnd_variant, Translocation)
        
        # Test basic properties
        self.assertEqual(bnd_variant.chrom, '20')
        self.assertEqual(bnd_variant.ref, 'G')
        self.assertEqual(bnd_variant.alt, 'G[2:321682[')
        
        # Check the destination properties only if it's a Translocation
        if isinstance(bnd_variant, Translocation):
            self.assertEqual(bnd_variant.dest_chrom, '2')
            self.assertEqual(bnd_variant.dest_pos, 321682)
    
    def test_record_properties(self):
        """Test VcfRecord properties that are used by variants."""
        # Test genotypes
        genotypes = self.snp_record.get_genotypes()
        self.assertEqual(len(genotypes), 2)
        self.assertEqual(genotypes[0].allele_indices, [0, 0])
        self.assertTrue(genotypes[0].phased)
        self.assertEqual(genotypes[1].allele_indices, [1, 0])
        self.assertTrue(genotypes[1].phased)
        
        # Test sample names
        self.assertEqual(self.snp_record.get_sample_names(), ['NA00001', 'NA00002'])
    
    def test_variant_string_representation(self):
        """Test the string representation of variants."""
        # Get variants using the iterator
        snp_variant = next(iter(self.snp_record))
        ins_variant = next(iter(self.ins_record))
        del_variant = next(iter(self.del_record))
        
        # Test string representation
        self.assertEqual(str(snp_variant), "SNP(rs6054257, 20:14369-14370, G>A)")
        self.assertEqual(str(ins_variant), "Insertion(variant_20_17329, 20:17329-17330, T>TA)")
        self.assertEqual(str(del_variant), "Deletion(rs6040355, 20:1110695-1110697, AT>A)")


if __name__ == '__main__':
    unittest.main()