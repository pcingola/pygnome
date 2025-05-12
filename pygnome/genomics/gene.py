"""Gene class for genomic annotations."""

from .biotype import Biotype
from .genomic_feature import GenomicFeature
from .transcript import Transcript


class Gene(GenomicFeature):
    """A gene, which may have multiple transcripts."""
    name: str | None = None
    biotype: Biotype | None = None
    transcripts: list[Transcript] = []
    
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