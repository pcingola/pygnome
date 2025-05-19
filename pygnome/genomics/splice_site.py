"""Splice site class for genomic annotations."""

from dataclasses import dataclass
from enum import Enum

from .genomic_feature import GenomicFeature


class SpliceSiteType(str, Enum):
    """Type of splice site."""
    DONOR = "donor"      # 5' splice site (GT)
    ACCEPTOR = "acceptor"  # 3' splice site (AG)


@dataclass
class SpliceSite(GenomicFeature):
    """A splice site (donor or acceptor) in a transcript."""
    site_type: SpliceSiteType