"""
Enumerations for VCF ANN field parsing.

This module defines enumerations used in parsing the ANN field in VCF files.
"""
from enum import Enum, auto


class AnnotationImpact(str, Enum):
    """
    Enumeration of putative impact categories for variant annotations.
    
    These impact categories represent the severity of the variant effect.
    """
    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    MODIFIER = "MODIFIER"


class FeatureType(str, Enum):
    """
    Enumeration of feature types for variant annotations.
    
    These types represent the genomic features affected by variants.
    """
    TRANSCRIPT = "transcript"
    MOTIF = "motif"
    MIRNA = "miRNA"
    GENE = "gene"
    INTERGENIC = "intergenic"
    EXON = "exon"
    INTRON = "intron"
    UPSTREAM = "upstream"
    DOWNSTREAM = "downstream"
    SPLICE_SITE = "splice_site"
    SPLICE_REGION = "splice_region"
    CHIP_SEQ_PEAK = "ChipSeq_peak"
    HISTONE_MARK = "histone_mark"
    CUSTOM = "custom"


class BiotypeCoding(str, Enum):
    """
    Enumeration of transcript biotype coding status.
    
    These values represent whether a transcript is coding or non-coding.
    """
    CODING = "Coding"
    NONCODING = "Noncoding"


class ErrorWarningType(str, Enum):
    """
    Enumeration of error and warning types in variant annotations.
    
    These codes represent various issues that can affect annotation accuracy.
    """
    # Error codes
    ERROR_CHROMOSOME_NOT_FOUND = "E1"
    ERROR_OUT_OF_CHROMOSOME_RANGE = "E2"
    
    # Warning codes
    WARNING_REF_DOES_NOT_MATCH_GENOME = "W1"
    WARNING_SEQUENCE_NOT_AVAILABLE = "W2"
    WARNING_TRANSCRIPT_INCOMPLETE = "W3"
    WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS = "W4"
    WARNING_TRANSCRIPT_NO_START_CODON = "W5"
    WARNING_TRANSCRIPT_NO_STOP_CODON = "W6"
    
    # Info codes
    INFO_REALIGN_3_PRIME = "I1"
    INFO_COMPOUND_ANNOTATION = "I2"
    INFO_NON_REFERENCE_ANNOTATION = "I3"
    
    @classmethod
    def from_code(cls, code: str) -> "ErrorWarningType | None":
        """
        Get the enum value from a code string.
        
        Args:
            code: The error/warning code (e.g., "E1", "W1", "I1")
            
        Returns:
            The corresponding ErrorWarningType enum value, or None if not found
        """
        for member in cls:
            if member.value == code:
                return member
        return None