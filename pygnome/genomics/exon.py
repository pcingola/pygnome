"""Exon class for genomic annotations."""

from dataclasses import dataclass

from .genomic_feature import GenomicFeature
from .phase import Phase


@dataclass
class Exon(GenomicFeature):
    """An exon within a transcript."""
    phase: Phase | None = None  # Frame phase for CDS; usually None for non-CDS exons