"""Splice site class for genomic annotations."""

from .genomic_feature import GenomicFeature
from .splice_site_type import SpliceSiteType


class SpliceSite(GenomicFeature):
    """A splice site (donor or acceptor) in a transcript."""
    site_type: SpliceSiteType