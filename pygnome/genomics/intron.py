"""Intron class for genomic annotations."""

from .genomic_feature import GenomicFeature
from .splice_site import SpliceSite


class Intron(GenomicFeature):
    """An intron within a transcript."""
    donor_site: SpliceSite | None = None
    acceptor_site: SpliceSite | None = None