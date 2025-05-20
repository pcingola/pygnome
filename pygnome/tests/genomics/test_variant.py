"""
Tests for the variant classes.
"""

import unittest
from pygnome.genomics.variant import (
    Variant, SNP, Insertion, Deletion, 
    Duplication, Inversion, Translocation, ComplexVariant
)
from pygnome.genomics.strand import Strand


class TestVariant(unittest.TestCase):
    """Test cases for the Variant class and its subclasses."""
    
    def test_base_variant(self):
        """Test creating a base Variant."""
        variant = Variant(
            id="test_variant",
            chrom="chr1",
            start=1000,
            end=1001,
            strand=Strand.POSITIVE,
            ref="A",
            alt="T"
        )
        
        self.assertEqual(variant.id, "test_variant")
        self.assertEqual(variant.chrom, "chr1")
        self.assertEqual(variant.start, 1000)
        self.assertEqual(variant.end, 1001)
        self.assertEqual(variant.strand, Strand.POSITIVE)
        self.assertEqual(variant.ref, "A")
        self.assertEqual(variant.alt, "T")
        self.assertEqual(variant.length, 1)
    
    
    def test_snp(self):
        """Test creating an SNP."""
        snp = SNP(
            id="rs123",
            chrom="chr1",
            start=1000,
            end=1001,
            strand=Strand.POSITIVE,
            ref="A",
            alt="T"
        )
        
        self.assertEqual(snp.ref, "A")
        self.assertEqual(snp.alt, "T")
        self.assertEqual(snp.length, 1)
    
    def test_invalid_snp(self):
        """Test that invalid SNPs raise validation errors."""
        # Reference allele too long
        with self.assertRaises(ValueError):
            SNP(
                id="invalid_ref",
                chrom="chr1",
                start=1000,
                end=1002,
                strand=Strand.POSITIVE,
                ref="AT",
                alt="T"
            )
        
        # Alternate allele too long
        with self.assertRaises(ValueError):
            SNP(
                id="invalid_alt",
                chrom="chr1",
                start=1000,
                end=1001,
                strand=Strand.POSITIVE,
                ref="A",
                alt="TG"
            )
    
    def test_insertion(self):
        """Test creating an Insertion."""
        insertion = Insertion(
            id="ins1",
            chrom="chr1",
            start=1000,
            end=1001,
            strand=Strand.POSITIVE,
            ref="A",
            alt="ATG"
        )
        
        self.assertEqual(insertion.ref, "A")
        self.assertEqual(insertion.alt, "ATG")
        self.assertEqual(insertion.inserted_sequence, "TG")
    
    def test_invalid_insertion(self):
        """Test that invalid Insertions raise validation errors."""
        # Alternate allele not longer than reference
        with self.assertRaises(ValueError):
            Insertion(
                id="invalid_ins",
                chrom="chr1",
                start=1000,
                end=1001,
                strand=Strand.POSITIVE,
                ref="A",
                alt="A"
            )
    
    def test_deletion(self):
        """Test creating a Deletion."""
        deletion = Deletion(
            id="del1",
            chrom="chr1",
            start=1000,
            end=1003,
            strand=Strand.POSITIVE,
            ref="ATG",
            alt="A"
        )
        
        self.assertEqual(deletion.ref, "ATG")
        self.assertEqual(deletion.alt, "A")
        self.assertEqual(deletion.deleted_sequence, "TG")
    
    def test_invalid_deletion(self):
        """Test that invalid Deletions raise validation errors."""
        # Alternate allele not shorter than reference
        with self.assertRaises(ValueError):
            Deletion(
                id="invalid_del",
                chrom="chr1",
                start=1000,
                end=1003,
                strand=Strand.POSITIVE,
                ref="ATG",
                alt="ATG"
            )
    
    def test_duplication(self):
        """Test creating a Duplication."""
        duplication = Duplication(
            id="dup1",
            chrom="chr1",
            start=1000,
            end=1003,
            strand=Strand.POSITIVE,
            ref="ATG",
            alt="<DUP>",
            dup_length=3
        )
        
        self.assertEqual(duplication.ref, "ATG")
        self.assertEqual(duplication.alt, "<DUP>")
        self.assertEqual(duplication.dup_length, 3)
        self.assertEqual(duplication.duplicated_sequence, "ATG")
    
    def test_inversion(self):
        """Test creating an Inversion."""
        inversion = Inversion(
            id="inv1",
            chrom="chr1",
            start=1000,
            end=1003,
            strand=Strand.POSITIVE,
            ref="ATG",
            alt="<INV>",
            inv_length=3
        )
        
        self.assertEqual(inversion.ref, "ATG")
        self.assertEqual(inversion.alt, "<INV>")
        self.assertEqual(inversion.inv_length, 3)
        self.assertEqual(inversion.inverted_sequence, "GTA")  # ATG reversed
    
    def test_translocation(self):
        """Test creating a Translocation."""
        translocation = Translocation(
            id="tra1",
            chrom="chr1",
            start=1000,
            end=1003,
            strand=Strand.POSITIVE,
            ref="ATG",
            alt="<TRA>",
            dest_chrom="chr2",
            dest_pos=2000
        )
        
        self.assertEqual(translocation.ref, "ATG")
        self.assertEqual(translocation.alt, "<TRA>")
        self.assertEqual(translocation.dest_chrom, "chr2")
        self.assertEqual(translocation.dest_pos, 2000)
    
    def test_complex_variant(self):
        """Test creating a ComplexVariant."""
        complex_var = ComplexVariant(
            id="complex1",
            chrom="chr1",
            start=1000,
            end=1010,
            strand=Strand.POSITIVE,
            ref="ATGCATGCAT",
            alt="ACGTACGTAC",
            description="Complex rearrangement"
        )
        
        self.assertEqual(complex_var.ref, "ATGCATGCAT")
        self.assertEqual(complex_var.alt, "ACGTACGTAC")
        self.assertEqual(complex_var.description, "Complex rearrangement")
    
    def test_variant_equality(self):
        """Test variant equality comparison."""
        # Create two identical variants
        variant1 = Variant(
            id="test_variant1",
            chrom="chr1",
            start=1000,
            end=1001,
            strand=Strand.POSITIVE,
            ref="A",
            alt="T"
        )
        
        variant2 = Variant(
            id="test_variant2",  # Different ID shouldn't affect equality
            chrom="chr1",
            start=1000,
            end=1001,
            strand=Strand.POSITIVE,
            ref="A",
            alt="T"
        )
        
        # Test equality
        self.assertEqual(variant1, variant2)
        
        # Create a variant with different case for ref and alt
        variant3 = Variant(
            id="test_variant3",
            chrom="chr1",
            start=1000,
            end=1001,
            strand=Strand.POSITIVE,
            ref="a",  # lowercase
            alt="t"   # lowercase
        )
        
        # Test equality with different case (should be equal because of case normalization in __post_init__)
        self.assertEqual(variant1, variant3)
        
        # Verify that the ref and alt were converted to uppercase
        self.assertEqual(variant3.ref, "A")
        self.assertEqual(variant3.alt, "T")
        
        # Create a variant with different attributes
        variant4 = Variant(
            id="test_variant4",
            chrom="chr1",
            start=1000,
            end=1001,
            strand=Strand.POSITIVE,
            ref="A",
            alt="G"  # Different alt
        )
        
        # Test inequality
        self.assertNotEqual(variant1, variant4)
    
    def test_variant_attribute_setting(self):
        """Test that setting ref and alt after initialization converts them to uppercase."""
        # Create a variant
        variant = Variant(
            id="test_variant",
            chrom="chr1",
            start=1000,
            end=1001,
            strand=Strand.POSITIVE,
            ref="A",
            alt="T"
        )
        
        # Set ref and alt to lowercase after initialization
        variant.ref = "g"
        variant.alt = "c"
        
        # Verify that they were converted to uppercase
        self.assertEqual(variant.ref, "G")
        self.assertEqual(variant.alt, "C")


if __name__ == "__main__":
    unittest.main()