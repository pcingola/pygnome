"""Gene class for genomic annotations."""

from typing import Any
from dataclasses import dataclass, field

from .biotype import Biotype
from .genomic_feature import GenomicFeature
from .transcript import Transcript


@dataclass
class Gene(GenomicFeature):
    """A gene, which may have multiple transcripts."""
    name: str | None = None
    biotype: Biotype | None = None
    transcripts: list[Transcript] = field(default_factory=list)
    chromosome: Any = None  # Reference to Chromosome
    
    def __post_init__(self):
        """Initialize the gene after creation."""
        super().__post_init__()
        # Initialize transcripts with references to this gene
        for transcript in self.transcripts:
            transcript.gene = self
    
    def add_transcript(self, transcript: Transcript) -> None:
        """Add a transcript to the gene."""
        self.transcripts.append(transcript)
        transcript.gene = self

    @property
    def canonical_transcript(self) -> Transcript | None:
        """Return the canonical transcript, if defined."""
        # In a real implementation, this would use some heuristic
        # like longest CDS or most exons, or a flag from the annotation
        if not self.transcripts:
            return None
        return self.transcripts[0]
    
    @property
    def is_coding(self) -> bool:
        """Return True if any transcript has coding sequence."""
        return any(transcript.is_coding for transcript in self.transcripts)
    
    def __iter__(self):
        """Iterate over transcripts sorted by start position."""
        return iter(sorted(self.transcripts, key=lambda x: x.start))

# No longer needed with dataclasses
# def update_forward_refs():
#     from .chromosome import Chromosome
#     Gene.model_rebuild()