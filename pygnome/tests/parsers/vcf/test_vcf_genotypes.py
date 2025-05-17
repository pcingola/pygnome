"""
Tests for the VcfGenotypes class.
"""
import unittest

from pygnome.parsers.vcf.vcf_header import VcfHeader
from pygnome.parsers.vcf.vcf_genotypes import VcfGenotypes, Genotype


class TestVcfGenotypes(unittest.TestCase):
    """Test cases for the VcfGenotypes class."""
    
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
        
        # Add sample names
        self.header.add_samples(['NA00001', 'NA00002', 'NA00003'])
        
        # Create fields for a VCF record
        self.fields = [
            '20',                   # CHROM
            '14370',                # POS
            'rs6054257',            # ID
            'G',                    # REF
            'A',                    # ALT
            '29',                   # QUAL
            'PASS',                 # FILTER
            'NS=3;DP=14;AF=0.5;DB', # INFO
            'GT:GQ:DP:HQ',          # FORMAT
            '0|0:48:1:51,51',       # Sample 1
            '1|0:48:8:51,51',       # Sample 2
            '1/1:43:5:.,.'          # Sample 3
        ]
        
        # Create a VcfGenotypes object with format and genotypes string
        format_and_genotypes_str = "GT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,."
        self.genotypes = VcfGenotypes(format_and_genotypes_str, self.header)
    
    def test_get_format(self):
        """Test getting the FORMAT field object."""
        format_obj = self.genotypes.get_format()
        self.assertIsNotNone(format_obj)
        keys = format_obj.get_keys() if format_obj else []
        self.assertEqual(keys, ["GT", "GQ", "DP", "HQ"])
    
    def test_get_format_keys(self):
        """Test getting the keys in the FORMAT field."""
        keys = self.genotypes.get_format_keys()
        self.assertEqual(keys, ["GT", "GQ", "DP", "HQ"])
    
    def test_has_format_key(self):
        """Test checking if a key is present in the FORMAT field."""
        self.assertTrue(self.genotypes.has_format_key("GT"))
        self.assertTrue(self.genotypes.has_format_key("GQ"))
        self.assertTrue(self.genotypes.has_format_key("DP"))
        self.assertTrue(self.genotypes.has_format_key("HQ"))
        self.assertFalse(self.genotypes.has_format_key("XX"))
    
    def test_get_value(self):
        """Test getting the value of a field for a specific sample."""
        # Check values for the first sample
        self.assertEqual(self.genotypes.get_value("GT", 0), "0|0")
        self.assertEqual(self.genotypes.get_value("GQ", 0), 48)
        self.assertEqual(self.genotypes.get_value("DP", 0), 1)
        
        # Parse the HQ field manually since it might not be parsed correctly in the test
        hq_value = self.genotypes.get_value("HQ", 0)
        self.assertIsNotNone(hq_value)
        self.assertEqual(hq_value, [51, 51])
        
        # Check values for the second sample
        self.assertEqual(self.genotypes.get_value("GT", 1), "1|0")
        self.assertEqual(self.genotypes.get_value("GQ", 1), 48)
        self.assertEqual(self.genotypes.get_value("DP", 1), 8)
        self.assertEqual(self.genotypes.get_value("HQ", 1), [51, 51])
        
        # Check values for the third sample
        self.assertEqual(self.genotypes.get_value("GT", 2), "1/1")
        self.assertEqual(self.genotypes.get_value("GQ", 2), 43)
        self.assertEqual(self.genotypes.get_value("DP", 2), 5)
        self.assertEqual(self.genotypes.get_value("HQ", 2), [None, None])  # .,.
    
    def test_set_value(self):
        """Test setting the value of a field for a specific sample."""
        # Set a value for the first sample
        self.genotypes.set_value("GQ", 50, 0)
        
        # Check that the value was set
        self.assertEqual(self.genotypes.get_value("GQ", 0), 50)
        
        # Set a value for the second sample
        self.genotypes.set_value("DP", 10, 1)
        
        # Check that the value was set
        self.assertEqual(self.genotypes.get_value("DP", 1), 10)
        
        # Set a value for the third sample
        self.genotypes.set_value("HQ", [60, 60], 2)
        
        # Check that the value was set
        self.assertEqual(self.genotypes.get_value("HQ", 2), [60, 60])
    
    def test_get_genotypes(self):
        """Test getting the genotypes for all samples."""
        genotypes = self.genotypes.get_genotypes()
        
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
    
    def test_set_genotype(self):
        """Test setting the genotype for a specific sample."""
        # Create a new genotype
        genotype = Genotype([0, 1], True)
        
        # Set the genotype for the first sample
        self.genotypes.set_genotype(genotype, 0)
        
        # Check that the genotype was set
        genotypes = self.genotypes.get_genotypes()
        self.assertEqual(genotypes[0].allele_indices, [0, 1])
        self.assertTrue(genotypes[0].phased)
        
        # Check that the GT field was updated
        self.assertEqual(self.genotypes.get_value("GT", 0), "0|1")
    
    def test_get_updated_format_and_genotypes(self):
        """Test getting updated format and genotypes string."""
        # Set a value for the first sample
        self.genotypes.set_value("GQ", 50, 0)
        
        # Get the updated format and genotypes string
        updated_str = self.genotypes.get_updated_format_and_genotypes()
        
        # Split it into fields
        fields = updated_str.split("\t")
        
        # Check the format field
        self.assertEqual(fields[0], "GT:GQ:DP:HQ")
        
        # Check the updated genotype field for the first sample
        self.assertEqual(fields[1], "0|0:50:1:51,51")
        
        # Check that the other genotype fields are unchanged
        self.assertEqual(fields[2], "1|0:48:8:51,51")
        self.assertEqual(fields[3], "1/1:43:5:.,.")
    
    def test_lazy_parsing(self):
        """Test lazy parsing of genotype data."""
        # Create a format and genotypes string with many samples
        format_and_genotypes_str = "GT:GQ:DP:HQ"
        for i in range(100):
            format_and_genotypes_str += f"\t0/0:{i}:1:51,51"
        
        # Add sample names
        header = VcfHeader()
        header.add_meta_line('##fileformat=VCFv4.5')
        header.add_meta_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        header.add_meta_line('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">')
        header.add_meta_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">')
        header.add_meta_line('##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">')
        
        # Add sample names
        sample_names = []
        for i in range(100):
            sample_names.append(f"SAMPLE{i}")
        header.add_samples(sample_names)
        
        # Create a VcfGenotypes object
        genotypes = VcfGenotypes(format_and_genotypes_str, header)

        # Verify that the fields haven't been parsed yet
        self.assertIsNone(genotypes._fields)

        # Get a value for a specific sample
        self.assertEqual(genotypes.get_value("GQ", 50), 50)  # 0-based index, sample 50 has GQ=50
        
        # Check that only the requested sample was parsed
        self.assertEqual(len(genotypes._sample_data_cache), 1)
        
        # Get a value for another sample
        self.assertEqual(genotypes.get_value("GQ", 75), 75)  # 0-based index, sample 75 has GQ=75
        
        # Check that only the requested samples were parsed
        self.assertEqual(len(genotypes._sample_data_cache), 2)
        
        # Verify that the fields have been parsed (get_value already calls _ensure_fields_parsed)
        self.assertIsNotNone(genotypes._fields)
        if genotypes._fields is not None:  # Add a check to avoid the type error
            self.assertEqual(len(genotypes._fields), 101)  # FORMAT + 100 samples
        
    def test_lazy_string_conversion(self):
        """Test that string conversion is lazy and uses the cached raw string if no changes."""
        # Create a format and genotypes string
        format_and_genotypes_str = "GT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,."
        
        # Create a VcfGenotypes object
        genotypes = VcfGenotypes(format_and_genotypes_str, self.header)
        
        # Get the updated format and genotypes string without making any changes
        updated_str = genotypes.get_updated_format_and_genotypes()
        
        # Verify that the original string is returned
        self.assertEqual(updated_str, format_and_genotypes_str)
        
        # Verify that the fields haven't been parsed
        self.assertIsNone(genotypes._fields)
        
        # Make a change
        genotypes.set_value("GQ", 50, 0)
        
        # Get the updated format and genotypes string
        updated_str = genotypes.get_updated_format_and_genotypes()
        
        # Verify that the fields have been parsed
        self.assertIsNotNone(genotypes._fields)
        
        # Verify that the string has been updated
        self.assertNotEqual(updated_str, format_and_genotypes_str)
        
        # Verify that the raw string has been invalidated
        self.assertIsNone(genotypes._raw_str)
        
    def test_str_and_repr(self):
        """Test the __str__ and __repr__ methods."""
        # Create a format and genotypes string
        format_and_genotypes_str = "GT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,."
        
        # Create a VcfGenotypes object
        genotypes = VcfGenotypes(format_and_genotypes_str, self.header)
        
        # Check that __str__ returns the format and genotypes string
        self.assertEqual(str(genotypes), format_and_genotypes_str)
        
        # Check that __repr__ returns the format and genotypes string
        self.assertEqual(repr(genotypes), format_and_genotypes_str)


if __name__ == '__main__':
    unittest.main()