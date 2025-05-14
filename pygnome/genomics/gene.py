"""Gene class for genomic annotations."""

from typing import Any, Dict, List, Optional, Union

from pydantic import Field

from .biotype import Biotype
from .genomic_feature import GenomicFeature
from .transcript import Transcript


class Gene(GenomicFeature):
    """A gene, which may have multiple transcripts."""
    name: str | None = None
    biotype: Biotype | None = None
    transcripts: list[Transcript] = []
    chromosome: Any = None  # Reference to Chromosome
    
    def __init__(self, **data):
        super().__init__(**data)
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

# This will be called after all models are defined
# to resolve forward references
def update_forward_refs():
    from .chromosome import Chromosome
    Gene.model_rebuild()