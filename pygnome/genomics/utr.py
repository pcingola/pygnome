"""UTR (Untranslated Region) class for genomic annotations."""

from enum import Enum
from dataclasses import dataclass

from .genomic_feature import GenomicFeature


class UTRType(str, Enum):
    """Type of UTR (untranslated region)."""
    FIVE_PRIME = "5'"
    THREE_PRIME = "3'"


@dataclass
class UTR(GenomicFeature):
    """An untranslated region (UTR) of a transcript."""
    utr_type: UTRType