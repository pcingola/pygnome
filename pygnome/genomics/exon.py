"""Exon class for genomic annotations."""

from pydantic import Field

from .genomic_feature import GenomicFeature
from .phase import Phase


class Exon(GenomicFeature):
    """An exon within a transcript."""
    phase: Phase | None = Field(
        None, description="Frame phase for CDS; usually None for non-CDS exons"
    )