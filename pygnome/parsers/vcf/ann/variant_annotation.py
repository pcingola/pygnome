"""
VariantAnnotation class for representing parsed ANN field data.

This module defines the VariantAnnotation class, which represents a single
annotation entry from the ANN field in a VCF record.
"""
from dataclasses import dataclass, field

from .enums import AnnotationImpact, FeatureType, ErrorWarningType, BiotypeCoding


@dataclass
class VariantAnnotation:
    """
    Represents a single annotation entry from the ANN field in a VCF record.
    
    This class encapsulates all the information contained in one annotation
    entry, following the VCF annotation format specification.
    """
    # Required fields
    allele: str
    annotation: str
    putative_impact: AnnotationImpact
    
    # Optional fields (may be None)
    gene_name: str | None = None
    gene_id: str | None = None
    feature_type: FeatureType | None = None
    feature_id: str | None = None
    transcript_biotype: BiotypeCoding | None = None
    
    # Positional information
    rank: int | None = None
    total: int | None = None
    
    # HGVS notation
    hgvs_c: str | None = None
    hgvs_p: str | None = None
    
    # Position information
    cdna_pos: int | None = None
    cdna_length: int | None = None
    cds_pos: int | None = None
    cds_length: int | None = None
    protein_pos: int | None = None
    protein_length: int | None = None
    
    # Distance to feature
    distance: int | None = None
    
    # Errors, warnings, or information messages
    messages: list[ErrorWarningType] | None = None
    
    def add_message(self, message: ErrorWarningType) -> None:
        """Add a message to the list of messages."""
        if self.messages is None:
            self.messages = []
        self.messages.append(message)

    def __str__(self) -> str:
        """Return a string representation of the annotation."""
        parts = [
            f"Allele: {self.allele}",
            f"Effect: {self.annotation}",
            f"Impact: {self.putative_impact.value}"
        ]
        
        if self.gene_name:
            parts.append(f"Gene: {self.gene_name}")
        
        if self.feature_type and self.feature_id:
            parts.append(f"{self.feature_type.value}: {self.feature_id}")
        
        if self.hgvs_c:
            parts.append(f"HGVS.c: {self.hgvs_c}")
        
        if self.hgvs_p:
            parts.append(f"HGVS.p: {self.hgvs_p}")
        
        if self.messages:
            msg_str = ", ".join(m.value for m in self.messages)
            parts.append(f"Messages: {msg_str}")
        
        return " | ".join(parts)