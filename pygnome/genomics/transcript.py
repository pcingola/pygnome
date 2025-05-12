"""Transcript class for genomic annotations."""

from .biotype import Biotype
from .cds import CDS
from .exon import Exon
from .genomic_feature import GenomicFeature
from .intron import Intron
from .splice_site import SpliceSite
from .utr import UTR


class Transcript(GenomicFeature):
    """A transcript of a gene, composed of exons, introns, UTRs, and CDS segments."""
    gene_id: str
    biotype: Biotype | None = None
    exons: list[Exon] = []
    cds_list: list[CDS] = []
    utrs: list[UTR] = []
    introns: list[Intron] = []
    splice_sites: list[SpliceSite] = []
    
    @property
    def coding_length(self) -> int:
        """Return the total length of coding sequence."""
        return sum(cds.length for cds in self.cds_list)
    
    @property
    def exonic_length(self) -> int:
        """Return the total length of exons."""
        return sum(exon.length for exon in self.exons)
    
    @property
    def is_coding(self) -> bool:
        """Return True if the transcript has coding sequence."""
        return len(self.cds_list) > 0
    
    @property
    def five_prime_utrs(self) -> list[UTR]:
        """Return all 5' UTRs."""
        from .utr_type import UTRType
        return [utr for utr in self.utrs if utr.utr_type == UTRType.FIVE_PRIME]
    
    @property
    def three_prime_utrs(self) -> list[UTR]:
        """Return all 3' UTRs."""
        from .utr_type import UTRType
        return [utr for utr in self.utrs if utr.utr_type == UTRType.THREE_PRIME]
    
    def __iter__(self):
        """Iterate over exons sorted by start position."""
        return iter(sorted(self.exons, key=lambda x: x.start))