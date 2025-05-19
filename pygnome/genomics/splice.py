"""Splice site and region classes for genomic annotations."""

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


@dataclass
class SpliceRegion(GenomicFeature):
    """
    A splice region in a transcript.
    
    Splice regions are the areas surrounding splice sites that may contain
    regulatory elements important for splicing. They typically extend a few
    bases into the exon and a few bases into the intron.
    """
    site_type: SpliceSiteType  # Using the same enum as SpliceSite
    exonic: bool  # Whether this region is in an exon (True) or intron (False)