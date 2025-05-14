"""Chromosome class for genomic annotations."""

from typing import Any, Dict, List, Optional, Union

from ..sequences.dna_string import DnaString
from .gene import Gene


class Chromosome:
    """A chromosome containing genes and optionally a DNA sequence."""
    
    def __init__(self, name: str, length: int | None = None,
                 sequence: DnaString | None = None):
        """Initialize a chromosome."""
        self.name = name
        self._length = length
        self.sequence = sequence
        self.genes: dict[str, Gene] = {}
        self.genome: Any = None  # Reference to Genome
        
    @property
    def length(self) -> int | None:
        """Return the chromosome length."""
        if self._length is not None:
            return self._length
        if self.sequence is not None:
            return len(self.sequence)
        return None
    
    def add_gene(self, gene: Gene) -> None:
        """Add a gene to the chromosome."""
        if gene.chrom != self.name:
            raise ValueError(f"Gene {gene.id} has chromosome {gene.chrom}, "
                            f"expected {self.name}")
        self.genes[gene.id] = gene
        
        # Set the chromosome reference in the gene
        gene.chromosome = self

    def add_genes(self, genes: list[Gene]) -> None:
        """Add multiple genes to the chromosome."""
        for gene in genes:
            self.add_gene(gene)

    def get_gene(self, gene_id: str) -> Gene | None:
        """Get a gene by ID."""
        return self.genes.get(gene_id)
    
    def get_genes_in_region(self, start: int, end: int) -> list[Gene]:
        """
        Get all genes that overlap a genomic region.
        
        Args:
            start: Start position (0-based, inclusive)
            end: End position (0-based, exclusive)
            
        Returns:
            List of genes that overlap the region
        """
        return [gene for gene in self.genes.values() 
                if gene.start <= end and gene.end >= start]
    
    def __str__(self) -> str:
        """Return a string representation of the chromosome."""
        length_str = f", length={self.length}" if self.length is not None else ""
        genes_str = f", {len(self.genes)} genes" if self.genes else ""
        return f"Chromosome({self.name}{length_str}{genes_str})"
    
    def __iter__(self):
        """Iterate over genes sorted by start position."""
        return iter(sorted(self.genes.values(), key=lambda x: x.start))

# This will be called after all models are defined
# to resolve forward references
def update_forward_refs():
    from .genome import Genome
    Chromosome.model_rebuild()