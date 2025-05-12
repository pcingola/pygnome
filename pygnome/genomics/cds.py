"""CDS (Coding Sequence) class for genomic annotations."""

from pydantic import Field

from .genomic_feature import GenomicFeature
from .phase import Phase


class CDS(GenomicFeature):
    """A coding sequence segment within a transcript."""
    phase: Phase = Field(..., description="Frame phase")