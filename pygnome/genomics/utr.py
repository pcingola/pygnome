"""UTR (Untranslated Region) class for genomic annotations."""

from dataclasses import dataclass

from .genomic_feature import GenomicFeature
from .utr_type import UTRType


@dataclass
class UTR(GenomicFeature):
    """An untranslated region (UTR) of a transcript."""
    utr_type: UTRType