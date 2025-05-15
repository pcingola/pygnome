"""
Variant classes for representing genomic variants.

This module defines a hierarchy of classes for representing genomic variants,
building on top of the GenomicFeature class.
"""

from typing import Optional, Any
from pydantic import Field, validator

from .genomic_feature import GenomicFeature
from .strand import Strand


class Variant(GenomicFeature):
    """
    Base class for all genomic variants.
    
    Extends GenomicFeature to represent a genomic variant with reference and
    alternate alleles.
    """
    ref: str = Field(..., description="Reference allele")
    alt: list[str] = Field(..., description="Alternate alleles")
    quality: Optional[float] = Field(None, description="Quality score")
    filters: list[str] = Field(default_factory=list, description="Filter flags")
    info: dict[str, Any] = Field(default_factory=dict, description="INFO fields")
    
    @validator('alt')
    def alt_must_not_be_empty(cls, v):
        """Validate that there is at least one alternate allele."""
        if not v:
            raise ValueError("At least one alternate allele is required")
        return v
    
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


class SNP(Variant):
    """
    Single Nucleotide Polymorphism (SNP) variant.
    
    A variant where a single nucleotide is changed.
    """
    @validator('ref')
    def ref_must_be_single_base(cls, v):
        """Validate that the reference allele is a single base."""
        if len(v) != 1:
            raise ValueError("Reference allele must be a single base for SNP")
        return v
    
    @validator('alt')
    def alt_must_be_single_bases(cls, v):
        """Validate that all alternate alleles are single bases."""
        for a in v:
            if len(a) != 1:
                raise ValueError("Alternate alleles must be single bases for SNP")
        return v


class Insertion(Variant):
    """
    Insertion variant.
    
    A variant where one or more nucleotides are inserted.
    """
    @validator('alt')
    def alt_must_be_longer_than_ref(cls, v, values):
        """Validate that all alternate alleles are longer than the reference."""
        if 'ref' in values:
            ref = values['ref']
            for a in v:
                if len(a) <= len(ref):
                    raise ValueError("Alternate alleles must be longer than reference for insertion")
        return v
    
    @property
    def inserted_sequence(self) -> list[str]:
        """Get the inserted sequence for each alternate allele."""
        return [a[len(self.ref):] for a in self.alt]


class Deletion(Variant):
    """
    Deletion variant.
    
    A variant where one or more nucleotides are deleted.
    """
    @validator('alt')
    def alt_must_be_shorter_than_ref(cls, v, values):
        """Validate that all alternate alleles are shorter than the reference."""
        if 'ref' in values:
            ref = values['ref']
            for a in v:
                if len(a) >= len(ref):
                    raise ValueError("Alternate alleles must be shorter than reference for deletion")
        return v
    
    @property
    def deleted_sequence(self) -> str:
        """Get the deleted sequence."""
        return self.ref[len(self.alt[0]):]


class Duplication(Variant):
    """
    Duplication variant.
    
    A variant where a segment of DNA is duplicated.
    """
    dup_length: int = Field(..., description="Length of the duplicated segment")
    
    @property
    def duplicated_sequence(self) -> str:
        """Get the duplicated sequence."""
        return self.ref[-self.dup_length:] if self.dup_length <= len(self.ref) else self.ref


class Inversion(Variant):
    """
    Inversion variant.
    
    A variant where a segment of DNA is reversed.
    """
    inv_length: int = Field(..., description="Length of the inverted segment")
    
    @property
    def inverted_sequence(self) -> str:
        """Get the inverted sequence."""
        return self.ref[:self.inv_length][::-1] if self.inv_length <= len(self.ref) else self.ref[::-1]


class Translocation(Variant):
    """
    Translocation variant.
    
    A variant where a segment of DNA is moved to a different location.
    """
    dest_chrom: str = Field(..., description="Destination chromosome")
    dest_pos: int = Field(..., description="Destination position")
    
    @validator('dest_pos')
    def dest_pos_must_be_non_negative(cls, v):
        """Validate that the destination position is non-negative."""
        if v < 0:
            raise ValueError("Destination position must be non-negative")
        return v


class ComplexVariant(Variant):
    """
    Complex variant.
    
    A variant that involves multiple types of changes and cannot be classified
    as a simple SNP, insertion, deletion, etc.
    """
    description: str = Field(..., description="Description of the complex variant")