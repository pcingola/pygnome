"""
Tests for the property-based implementation of the Variant class.

These tests specifically verify that the ref and alt properties
correctly enforce uppercase values.
"""

import unittest
from pygnome.genomics.variant import (
    Variant, SNP, Insertion, Deletion, 
    Duplication, Inversion, Translocation, ComplexVariant
)
from pygnome.genomics.strand import Strand


class TestVariantProperties(unittest.TestCase):
    """Test cases for the property-based implementation of the Variant class."""
    
    def test_initialization_uppercase_conversion(self):
        """Test that ref and alt are converted to uppercase during initialization."""
        variant = Variant(
            id="test_variant",
            chrom="chr1",
            start=1000,
            end=1001,
            strand=Strand.POSITIVE,
            ref="a",  # lowercase
            alt="t"   # lowercase
        )
        
        # Verify that ref and alt were converted to uppercase
        self.assertEqual(variant.ref, "A")
        self.assertEqual(variant.alt, "T")
    
    def test_property_setter_uppercase_conversion(self):
        """Test that setting ref and alt after initialization converts them to uppercase."""
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
    
    def test_empty_values(self):
        """Test handling of empty values for ref and alt."""
        # Empty string should remain empty
        variant = Variant(
            id="test_empty",
            chrom="chr1",
            start=1000,
            end=1001,
            strand=Strand.POSITIVE,
            ref="",
            alt="A"  # Need a valid alt to avoid validation error
        )
        self.assertEqual(variant.ref, "")
        
        # Setting to empty string after initialization
        variant.ref = ""
        self.assertEqual(variant.ref, "")
    
    def test_derived_classes(self):
        """Test that properties work correctly with derived classes."""
        # Test SNP
        snp = SNP(
            id="rs123",
            chrom="chr1",
            start=1000,
            end=1001,
            strand=Strand.POSITIVE,
            ref="a",  # lowercase
            alt="t"   # lowercase
        )
        self.assertEqual(snp.ref, "A")
        self.assertEqual(snp.alt, "T")
        
        # Test Insertion
        insertion = Insertion(
            id="ins1",
            chrom="chr1",
            start=1000,
            end=1001,
            strand=Strand.POSITIVE,
            ref="a",  # lowercase
            alt="atg"  # lowercase
        )
        self.assertEqual(insertion.ref, "A")
        self.assertEqual(insertion.alt, "ATG")
        self.assertEqual(insertion.inserted_sequence, "TG")
        
        # Test Deletion
        deletion = Deletion(
            id="del1",
            chrom="chr1",
            start=1000,
            end=1003,
            strand=Strand.POSITIVE,
            ref="atg",  # lowercase
            alt="a"     # lowercase
        )
        self.assertEqual(deletion.ref, "ATG")
        self.assertEqual(deletion.alt, "A")
        self.assertEqual(deletion.deleted_sequence, "TG")
    
    def test_mixed_case(self):
        """Test handling of mixed case values for ref and alt."""
        variant = Variant(
            id="test_mixed",
            chrom="chr1",
            start=1000,
            end=1001,
            strand=Strand.POSITIVE,
            ref="AtG",  # mixed case
            alt="cGa"   # mixed case
        )
        
        # Verify that ref and alt were converted to uppercase
        self.assertEqual(variant.ref, "ATG")
        self.assertEqual(variant.alt, "CGA")
        
        # Set ref and alt to mixed case after initialization
        variant.ref = "gCt"
        variant.alt = "tAc"
        
        # Verify that they were converted to uppercase
        self.assertEqual(variant.ref, "GCT")
        self.assertEqual(variant.alt, "TAC")
    
    def test_equality_with_properties(self):
        """Test that equality comparison works correctly with properties."""
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
            id="test_variant2",
            chrom="chr1",
            start=1000,
            end=1001,
            strand=Strand.POSITIVE,
            ref="a",  # lowercase
            alt="t"   # lowercase
        )
        
        # Test equality with different case (should be equal because of case normalization)
        self.assertEqual(variant1, variant2)
        
        # Modify ref and alt of variant2
        variant2.ref = "g"
        variant2.alt = "c"
        
        # Test inequality after modification
        self.assertNotEqual(variant1, variant2)
    
    def test_property_persistence(self):
        """Test that property values persist correctly."""
        variant = Variant(
            id="test_persistence",
            chrom="chr1",
            start=1000,
            end=1001,
            strand=Strand.POSITIVE,
            ref="A",
            alt="T"
        )
        
        # Modify ref and alt
        variant.ref = "G"
        variant.alt = "C"
        
        # Verify that the changes persist
        self.assertEqual(variant.ref, "G")
        self.assertEqual(variant.alt, "C")
        
        # Modify again
        variant.ref = "atg"
        variant.alt = "cga"
        
        # Verify that the changes persist and are uppercase
        self.assertEqual(variant.ref, "ATG")
        self.assertEqual(variant.alt, "CGA")


if __name__ == "__main__":
    unittest.main()