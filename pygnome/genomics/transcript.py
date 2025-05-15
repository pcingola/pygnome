"""Transcript class for genomic annotations."""

from typing import Any
from dataclasses import dataclass, field

from .biotype import Biotype
from .cds import CDS
from .codon_table import CodonTable, CodonTableType
from .exon import Exon
from .genomic_feature import GenomicFeature
from .intron import Intron
from .splice_site import SpliceSite
from .strand import Strand
from .utr import UTR
from ..sequences.utils import reverse_complement


@dataclass
class Transcript(GenomicFeature):
    """A transcript of a gene, composed of exons, introns, UTRs, and CDS segments."""
    gene_id: str
    biotype: Biotype | None = None
    exons: list[Exon] = field(default_factory=list)
    cds_list: list[CDS] = field(default_factory=list)
    utrs: list[UTR] = field(default_factory=list)
    introns: list[Intron] = field(default_factory=list)
    splice_sites: list[SpliceSite] = field(default_factory=list)
    _cds: str | None = None
    _protein: str | None = None
    gene: Any = None  # Reference to Gene
    
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
    
    def cds(self, genome_sequence: dict[str, str]) -> str:
        """Get the coding DNA sequence for this transcript."""
        # Use cached value if available
        if self._cds is not None:
            return self._cds
        
        if not self.is_coding:
            self._cds = ""
            return ""
        
        # Get the chromosome sequence
        if self.chrom not in genome_sequence:
            raise ValueError(f"Chromosome {self.chrom} not found in genome sequence")
        
        chrom_seq = genome_sequence[self.chrom]
        
        # Sort CDS segments by position based on strand
        if self.strand == Strand.POSITIVE:
            # For positive strand, sort by increasing position
            sorted_cds = sorted(self.cds_list, key=lambda x: x.start)
        else:
            # For negative strand, sort by decreasing position
            sorted_cds = sorted(self.cds_list, key=lambda x: x.start, reverse=True)
        
        # Extract the coding sequence for each CDS segment
        coding_seq = ""
        for cds in sorted_cds:
            # Extract the sequence for this CDS segment
            segment_seq = chrom_seq[cds.start:cds.end]
            
            # Reverse complement if on negative strand
            if self.strand == Strand.NEGATIVE:
                segment_seq = reverse_complement(segment_seq)
                
            coding_seq += segment_seq
        
        # Cache the result
        self._cds = coding_seq
        
        return coding_seq
    
    def protein(self, genome_sequence: dict[str, str]) -> str:
        """Get the protein sequence for this transcript."""
        # Get the coding sequence
        coding_seq = self.cds(genome_sequence)
        
        if not coding_seq:
            return ""
        
        # Use the genome's codon table
        table_type = self.gene.chromosome.genome.codon_table_type
        
        # Create a codon table
        codon_table = CodonTable(table_type)
        
        # Translate the coding sequence
        protein_seq = codon_table.translate_sequence(coding_seq)
        
        return protein_seq

# No longer needed with dataclasses
# def update_forward_refs():
#     from .gene import Gene
#     Transcript.model_rebuild()