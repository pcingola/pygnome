"""
VCF Variant class for representing variants from VCF records.
"""
from enum import Enum, auto
from pygnome.parsers.vcf.vcf_record import VcfRecord, Genotype


class VariantType(str, Enum):
    """Enum representing different types of variants."""
    SNP = "SNP"
    INS = "INS"
    DEL = "DEL"
    SV = "SV"
    BND = "BND"
    OTHER = "OTHER"

class VcfVariant:
    """
    Represents a variant from a VCF record.
    
    This class provides a higher-level representation of a variant, with methods
    for determining the variant type and accessing alleles and genotypes.
    """
    def __init__(self, record: VcfRecord):
        self.record = record
    
    def get_chrom(self) -> str:
        """Get the chromosome name."""
        return self.record.get_chrom()
    
    def get_pos(self) -> int:
        """Get the 0-based position."""
        return self.record.get_pos()
    
    def get_end(self) -> int:
        """
        Get the end position of the variant (0-based, exclusive).
        
        For SNPs, this is pos + 1. For indels, it's pos + len(ref).
        For structural variants, it's determined by the END or SVLEN INFO field.
        """
        if self.is_structural_variant():
            # Check for END info field
            if self.record.has_info("END"):
                end = self.record.get_info("END")
                if end is not None:
                    # END is 1-based inclusive in VCF, convert to 0-based exclusive
                    return end
            
            # Check for SVLEN info field
            if self.record.has_info("SVLEN"):
                svlen = self.record.get_info("SVLEN")
                if svlen is not None:
                    # SVLEN is the length of the variant
                    if isinstance(svlen, list):
                        # Use the first SVLEN value
                        return self.get_pos() + abs(svlen[0])
                    else:
                        return self.get_pos() + abs(svlen)
        
        # Default: pos + len(ref)
        return self.get_pos() + len(self.record.get_ref())
    
    def get_id(self) -> str:
        """Get the ID field."""
        return self.record.get_id()
    
    def get_ref(self) -> str:
        """Get the reference allele."""
        return self.record.get_ref()
    
    def get_alt(self) -> list[str]:
        """Get the list of alternate alleles."""
        return self.record.get_alt()
    
    def get_alleles(self) -> list[str]:
        """Get all alleles (reference and alternates)."""
        return self.record.get_alleles()
    
    def get_genotypes(self) -> list[Genotype]:
        """Get the genotypes for all samples."""
        return self.record.get_genotypes()
    
    def get_sample_names(self) -> list[str]:
        """Get the names of the samples."""
        return self.record.get_sample_names()
    
    def is_snp(self) -> bool:
        """
        Check if this variant is a SNP (Single Nucleotide Polymorphism).
        
        Returns:
            True if the variant is a SNP, False otherwise
        """
        ref = self.get_ref()
        alt = self.get_alt()
        
        # SNP: single reference base and all alternate alleles are single bases
        return (len(ref) == 1 and 
                all(len(a) == 1 for a in alt) and 
                all(a in "ACGTN" for a in alt))
    
    def is_indel(self) -> bool:
        """
        Check if this variant is an indel (insertion or deletion).
        
        Returns:
            True if the variant is an indel, False otherwise
        """
        # Structural variants are not considered indels
        if self.is_structural_variant() or self.is_breakend():
            return False
            
        ref = self.get_ref()
        alt = self.get_alt()
        
        # Indel: reference and at least one alternate allele have different lengths
        return any(len(ref) != len(a) for a in alt)
    
    def is_insertion(self) -> bool:
        """
        Check if this variant is an insertion.
        
        Returns:
            True if the variant is an insertion, False otherwise
        """
        # Structural variants are not considered insertions
        if self.is_structural_variant() or self.is_breakend():
            return False
            
        ref = self.get_ref()
        alt = self.get_alt()
        
        # Insertion: at least one alternate allele is longer than the reference
        return any(len(a) > len(ref) for a in alt)
    
    def is_deletion(self) -> bool:
        """
        Check if this variant is a deletion.
        
        Returns:
            True if the variant is a deletion, False otherwise
        """
        # Structural variants are not considered deletions
        if self.is_structural_variant() or self.is_breakend():
            return False
            
        ref = self.get_ref()
        alt = self.get_alt()
        
        # Deletion: at least one alternate allele is shorter than the reference
        return any(len(a) < len(ref) for a in alt)
    
    def is_structural_variant(self) -> bool:
        """
        Check if this variant is a structural variant.
        
        Returns:
            True if the variant is a structural variant, False otherwise
        """
        alt = self.get_alt()
        
        # Structural variant: at least one alternate allele is symbolic (<...>)
        return any(a.startswith("<") and a.endswith(">") for a in alt)
    
    def is_breakend(self) -> bool:
        """
        Check if this variant is a breakend.
        
        Returns:
            True if the variant is a breakend, False otherwise
        """
        alt = self.get_alt()
        
        # Breakend: at least one alternate allele contains [ or ]
        return any("[" in a or "]" in a for a in alt)
    
    def is_multi_allelic(self) -> bool:
        """
        Check if this variant has multiple alternate alleles.
        
        Returns:
            True if the variant has multiple alternate alleles, False otherwise
        """
        return len(self.get_alt()) > 1
    
    def get_variant_type(self) -> VariantType:
        """
        Get the type of this variant.
        
        Returns:
            A string describing the variant type (SNP, indel, structural variant, etc.)
        """
        if self.is_snp():
            return VariantType.SNP
        elif self.is_insertion():
            return VariantType.INS
        elif self.is_deletion():
            return VariantType.DEL
        elif self.is_structural_variant():
            # Get the specific type from the symbolic allele
            alt = self.get_alt()
            for a in alt:
                if a.startswith("<") and a.endswith(">"):
                    # Try to map the symbolic allele to a variant type
                    sv_type = a[1:-1]  # Remove the < and >
                    return sv_type
            return VariantType.SV
        elif self.is_breakend():
            return VariantType.BND
        else:
            return VariantType.OTHER
    
    def __str__(self) -> str:
        """Return a string representation of the variant."""
        return f"{self.get_chrom()}:{self.get_pos() + 1}:{self.get_ref()}>{','.join(self.get_alt())}"