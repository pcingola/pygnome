"""
Factory class for creating Variant objects from VCF records.
"""
from typing import Iterator

from pygnome.genomics.strand import Strand
from pygnome.genomics.variant import (
    VariantType, Variant, SNP, Insertion, Deletion,
    Duplication, Inversion, Translocation, ComplexVariant
)


class VariantFactory:
    """
    Factory class for creating Variant objects from VCF records.
    
    This class encapsulates the logic for determining variant types and
    creating appropriate Variant objects from VCF record data.
    """
    
    @staticmethod
    def is_snp(ref: str, alt: str) -> bool:
        """
        Check if a variant is a SNP (Single Nucleotide Polymorphism).
        
        Args:
            ref: Reference allele
            alt: Alternate allele
            
        Returns:
            True if the variant is a SNP, False otherwise
        """
        # SNP: single reference base and alternate allele is a single base
        return (len(ref) == 1 and len(alt) == 1 and alt in "ACGTN")
    
    @staticmethod
    def is_indel(ref: str, alt: str) -> bool:
        """
        Check if a variant is an indel (insertion or deletion).
        
        Args:
            ref: Reference allele
            alt: Alternate allele
            
        Returns:
            True if the variant is an indel, False otherwise
        """
        # Structural variants are not considered indels
        if VariantFactory.is_structural_variant(alt) or VariantFactory.is_breakend(alt):
            return False
            
        # Indel: reference and alternate allele have different lengths
        return len(ref) != len(alt)
    
    @staticmethod
    def is_insertion(ref: str, alt: str) -> bool:
        """
        Check if a variant is an insertion.
        
        Args:
            ref: Reference allele
            alt: Alternate allele
            
        Returns:
            True if the variant is an insertion, False otherwise
        """
        # Structural variants are not considered insertions
        if VariantFactory.is_structural_variant(alt) or VariantFactory.is_breakend(alt):
            return False
            
        # Insertion: alternate allele is longer than the reference
        return len(alt) > len(ref)
    
    @staticmethod
    def is_deletion(ref: str, alt: str) -> bool:
        """
        Check if a variant is a deletion.
        
        Args:
            ref: Reference allele
            alt: Alternate allele
            
        Returns:
            True if the variant is a deletion, False otherwise
        """
        # Structural variants are not considered deletions
        if VariantFactory.is_structural_variant(alt) or VariantFactory.is_breakend(alt):
            return False
            
        # Deletion: alternate allele is shorter than the reference
        return len(alt) < len(ref)
    
    @staticmethod
    def is_structural_variant(alt: str) -> bool:
        """
        Check if a variant is a structural variant.
        
        Args:
            alt: Alternate allele
            
        Returns:
            True if the variant is a structural variant, False otherwise
        """
        # Structural variant: alternate allele is symbolic (<...>)
        return alt.startswith("<") and alt.endswith(">")
    
    @staticmethod
    def is_breakend(alt: str) -> bool:
        """
        Check if a variant is a breakend.
        
        Args:
            alt: Alternate allele
            
        Returns:
            True if the variant is a breakend, False otherwise
        """
        # Breakend: alternate allele contains [ or ]
        return "[" in alt or "]" in alt
    
    @staticmethod
    def get_variant_type(ref: str, alt: str) -> VariantType:
        """
        Get the type of a variant.
        
        Args:
            ref: Reference allele
            alt: Alternate allele
            
        Returns:
            A VariantType enum value describing the variant type
        """
        if VariantFactory.is_snp(ref, alt):
            return VariantType.SNP
        elif VariantFactory.is_insertion(ref, alt):
            return VariantType.INS
        elif VariantFactory.is_deletion(ref, alt):
            return VariantType.DEL
        elif VariantFactory.is_structural_variant(alt):
            # Get the specific type from the symbolic allele
            if alt.startswith("<") and alt.endswith(">"):
                # Return SV type but preserve the original type in the string representation
                return VariantType.SV
            return VariantType.SV
        elif VariantFactory.is_breakend(alt):
            return VariantType.BND
        else:
            return VariantType.OTHER
    
    @staticmethod
    def create_variants_from_record(record) -> Iterator[Variant]:
        """
        Create Variant objects from a VCF record.
        
        Args:
            record: A VcfRecord object
            
        Yields:
            Variant objects for each alternate allele in the record
        """
        for i in range(len(record.alts)):
            variant_id = record.id if record.id else f"variant_{record.chrom}_{record.start}"
            alt = record.alts[i]
            
            # Basic variant parameters
            variant_params = {
                "id": variant_id,
                "chrom": record.chrom,
                "start": record.start,
                "end": record.get_end(i),
                "strand": Strand.UNSTRANDED,
                "ref": record.ref,
                "alt": alt
            }
            
            if VariantFactory.is_snp(record.ref, alt):
                yield SNP(**variant_params)
            elif VariantFactory.is_insertion(record.ref, alt):
                yield Insertion(**variant_params)
            elif VariantFactory.is_deletion(record.ref, alt):
                yield Deletion(**variant_params)
            elif VariantFactory.is_structural_variant(alt):
                # Handle specific structural variant types
                sv_type = alt[1:-1] if alt.startswith("<") and alt.endswith(">") else ""
                
                if sv_type == "DEL":
                    # Structural deletion - use the Deletion class
                    yield Deletion(**variant_params)
                elif sv_type == "DUP" or sv_type == "DUP:TANDEM":
                    # Duplication
                    dup_length = abs(record.get_info("SVLEN")) if record.has_info("SVLEN") else record.get_end(i) - record.get_pos()
                    yield Duplication(**variant_params, dup_length=dup_length)
                elif sv_type == "INV":
                    # Inversion
                    inv_length = abs(record.get_info("SVLEN")) if record.has_info("SVLEN") else record.get_end(i) - record.get_pos()
                    yield Inversion(**variant_params, inv_length=inv_length)
                elif sv_type == "TRA" or sv_type == "CTX":
                    # Translocation - requires additional info
                    dest_chrom = record.get_info("CHR2") if record.has_info("CHR2") else record.get_chrom()
                    dest_pos = record.get_info("POS2") if record.has_info("POS2") else 0
                    yield Translocation(**variant_params, dest_chrom=dest_chrom, dest_pos=dest_pos)
                else:
                    # Other structural variants
                    yield ComplexVariant(**variant_params, description=f"Structural variant ({sv_type})")
            elif VariantFactory.is_breakend(alt):
                # Parse breakend notation to extract destination
                # Breakend format: X[chr:pos[ or X]chr:pos] or [chr:pos[X or ]chr:pos]X
                bnd_alt = alt
                dest_chrom = ""
                dest_pos = 0
                
                # Extract destination from breakend notation
                if "[" in bnd_alt:
                    parts = bnd_alt.split("[")
                    for part in parts:
                        if ":" in part:
                            dest_info = part.split(":")
                            if len(dest_info) == 2:
                                dest_chrom = dest_info[0]
                                dest_pos = int(dest_info[1]) if dest_info[1].isdigit() else 0
                                break
                elif "]" in bnd_alt:
                    parts = bnd_alt.split("]")
                    for part in parts:
                        if ":" in part:
                            dest_info = part.split(":")
                            if len(dest_info) == 2:
                                dest_chrom = dest_info[0]
                                dest_pos = int(dest_info[1]) if dest_info[1].isdigit() else 0
                                break
                
                if dest_chrom and dest_pos > 0:
                    yield Translocation(**variant_params, dest_chrom=dest_chrom, dest_pos=dest_pos)
                else:
                    yield ComplexVariant(**variant_params, description=f"Breakend variant")
            else:
                # Fallback for truly complex variants that don't fit other categories
                yield ComplexVariant(**variant_params, description=f"Complex variant ({VariantFactory.get_variant_type(record.ref, alt)})")