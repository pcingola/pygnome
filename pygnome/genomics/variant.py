"""
Variant classes for representing genomic variants.

This module defines a hierarchy of classes for representing genomic variants,
building on top of the GenomicFeature class.
"""

from typing import Any
from dataclasses import dataclass, field

from .genomic_feature import GenomicFeature


@dataclass
class Variant(GenomicFeature):
    """
    Base class for all genomic variants.
    
    Extends GenomicFeature to represent a genomic variant with reference and
    alternate alleles.
    """
    ref: str  # Reference allele
    alt: list[str]  # Alternate alleles
    quality: float | None = None  # Quality score
    filters: list[str] = field(default_factory=list)  # Filter flags
    info: dict[str, Any] = field(default_factory=dict)  # INFO fields
    
    def __post_init__(self):
        """Validate the variant after initialization."""
        super().__post_init__()
        
        # Validate that there is at least one alternate allele
        if not self.alt:
            raise ValueError("At least one alternate allele is required")
    
    @property
    def is_multi_allelic(self) -> bool:
        """Check if this variant has multiple alternate alleles."""
        return len(self.alt) > 1
    
    @property
    def alleles(self) -> list[str]:
        """Get all alleles (reference and alternates)."""
        return [self.ref] + self.alt
    
    def __str__(self) -> str:
        """Return a string representation of the variant."""
        return f"{self.__class__.__name__}({self.id}, {self.chrom}:{self.start}-{self.end}, {self.ref}>{','.join(self.alt)})"


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
        
        # Validate that all alternate alleles are single bases
        for a in self.alt:
            if len(a) != 1:
                raise ValueError("Alternate alleles must be single bases for SNP")


@dataclass
class Insertion(Variant):
    """
    Insertion variant.
    
    A variant where one or more nucleotides are inserted.
    """
    
    def __post_init__(self):
        """Validate the insertion after initialization."""
        super().__post_init__()
        
        # Validate that all alternate alleles are longer than the reference
        for a in self.alt:
            if len(a) <= len(self.ref):
                raise ValueError("Alternate alleles must be longer than reference for insertion")
    
    @property
    def inserted_sequence(self) -> list[str]:
        """Get the inserted sequence for each alternate allele."""
        return [a[len(self.ref):] for a in self.alt]


@dataclass
class Deletion(Variant):
    """
    Deletion variant.
    
    A variant where one or more nucleotides are deleted.
    """
    
    def __post_init__(self):
        """Validate the deletion after initialization."""
        super().__post_init__()
        
        # Validate that all alternate alleles are shorter than the reference
        for a in self.alt:
            if len(a) >= len(self.ref):
                raise ValueError("Alternate alleles must be shorter than reference for deletion")
    
    @property
    def deleted_sequence(self) -> str:
        """Get the deleted sequence."""
        return self.ref[len(self.alt[0]):]


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