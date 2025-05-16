"""
Variant classes for representing genomic variants.

This module defines a hierarchy of classes for representing genomic variants,
building on top of the GenomicFeature class.
"""
from typing import Any
from dataclasses import dataclass, field
from enum import Enum, auto

from .genomic_feature import GenomicFeature


class VariantType(str, Enum):
    """Enum representing different types of variants."""
    SNP = "SNP"
    INS = "INS"
    DEL = "DEL"
    SV = "SV"
    BND = "BND"
    OTHER = "OTHER"



@dataclass
class Variant(GenomicFeature):
    """
    Base class for all genomic variants.
    
    Extends GenomicFeature to represent a genomic variant with reference and
    alternate alleles.
    """
    ref: str  # Reference allele
    alt: str  # Alternate allele
    
    def __post_init__(self):
        """Validate the variant after initialization."""
        super().__post_init__()
        
        # Validate that the alternate allele is not empty
        if not self.alt:
            raise ValueError("Alternate allele is required")
    def __str__(self) -> str:
        """Return a string representation of the variant."""
        return f"{self.__class__.__name__}({self.id}, {self.chrom}:{self.start}-{self.end}, {self.ref}>{self.alt})"


@dataclass
class SNP(Variant):
    """
    Single Nucleotide Polymorphism (SNP) variant.
    
    A variant where a single nucleotide is changed.
    """
    
    def __post_init__(self):
        """Validate the SNP after initialization."""
        super().__post_init__()
        
        # Validate that the reference allele is a single base
        if len(self.ref) != 1:
            raise ValueError("Reference allele must be a single base for SNP")
        
        # Validate that the alternate allele is a single base
        if len(self.alt) != 1:
            raise ValueError("Alternate allele must be a single base for SNP")


@dataclass
class Insertion(Variant):
    """
    Insertion variant.
    
    A variant where one or more nucleotides are inserted.
    """
    
    def __post_init__(self):
        """Validate the insertion after initialization."""
        super().__post_init__()
        
        # Validate that the alternate allele is longer than the reference
        if len(self.alt) <= len(self.ref):
            raise ValueError("Alternate allele must be longer than reference for insertion")
    
    @property
    def inserted_sequence(self) -> str:
        """Get the inserted sequence."""
        return self.alt[len(self.ref):]


@dataclass
class Deletion(Variant):
    """
    Deletion variant.
    
    A variant where one or more nucleotides are deleted.
    """
    
    def __post_init__(self):
        """Validate the deletion after initialization."""
        super().__post_init__()
        
        # Validate that the alternate allele is shorter than the reference
        # Skip validation for symbolic alleles like <DEL>
        if not (self.alt.startswith("<") and self.alt.endswith(">")):
            if len(self.alt) >= len(self.ref):
                raise ValueError("Alternate allele must be shorter than reference for deletion")
    
    @property
    def deleted_sequence(self) -> str:
        """Get the deleted sequence."""
        # For symbolic alleles, return the reference sequence
        if self.alt.startswith("<") and self.alt.endswith(">"):
            return self.ref
        return self.ref[len(self.alt):]


@dataclass
class Duplication(Variant):
    """
    Duplication variant.
    
    A variant where a segment of DNA is duplicated.
    """
    dup_length: int = field(default=0)  # Length of the duplicated segment
    
    def __post_init__(self):
        """Validate the duplication after initialization."""
        super().__post_init__()
        
        # Ensure dup_length is provided
        if self.dup_length <= 0:
            raise ValueError("dup_length must be a positive integer")
    
    @property
    def duplicated_sequence(self) -> str:
        """Get the duplicated sequence."""
        return self.ref[-self.dup_length:] if self.dup_length <= len(self.ref) else self.ref


@dataclass
class Inversion(Variant):
    """
    Inversion variant.
    
    A variant where a segment of DNA is reversed.
    """
    inv_length: int = field(default=0)  # Length of the inverted segment
    
    def __post_init__(self):
        """Validate the inversion after initialization."""
        super().__post_init__()
        
        # Ensure inv_length is provided
        if self.inv_length <= 0:
            raise ValueError("inv_length must be a positive integer")
    
    @property
    def inverted_sequence(self) -> str:
        """Get the inverted sequence."""
        return self.ref[:self.inv_length][::-1] if self.inv_length <= len(self.ref) else self.ref[::-1]


@dataclass
class Translocation(Variant):
    """
    Translocation variant.
    
    A variant where a segment of DNA is moved to a different location.
    """
    dest_chrom: str = field(default="")  # Destination chromosome
    dest_pos: int = field(default=0)     # Destination position
    
    def __post_init__(self):
        """Validate the translocation after initialization."""
        super().__post_init__()
        
        # Ensure dest_chrom is provided
        if not self.dest_chrom:
            raise ValueError("dest_chrom must be provided")
        
        # Validate that the destination position is non-negative
        if self.dest_pos < 0:
            raise ValueError("Destination position must be non-negative")


@dataclass
class ComplexVariant(Variant):
    """
    Complex variant.
    
    A variant that involves multiple types of changes and cannot be classified
    as a simple SNP, insertion, deletion, etc.
    """
    description: str = field(default="")  # Description of the complex variant
    
    def __post_init__(self):
        """Validate the complex variant after initialization."""
        super().__post_init__()
        
        # Ensure description is provided
        if not self.description:
            raise ValueError("description must be provided")