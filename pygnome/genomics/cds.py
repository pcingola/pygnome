"""CDS (Coding Sequence) class for genomic annotations."""

from dataclasses import dataclass

from .genomic_feature import GenomicFeature
from .phase import Phase


@dataclass
class CDS(GenomicFeature):
    """A coding sequence segment within a transcript."""
    phase: Phase  # Frame phase